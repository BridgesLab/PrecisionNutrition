# robust_mr_helpers.R -----------------------------------------------------------
# Helpers for the overlap- and pleiotropy-robust MR arm (CAUSE + MR-APSS) of the
# cholesterol -> BMD analysis. Sourced by robust_mr_*.qmd.
#
# Design notes:
#   * Everything is column-spec driven via config_robust_mr.yml, so a change in an
#     upstream file's header is a config edit, not a code edit.
#   * QC follows the MR-APSS published criteria (HapMap3, MAF>0.05, unambiguous
#     alleles, INFO>0.9, no MHC, chi2 < max(N/1000, 80)) so that CAUSE and MR-APSS
#     see the *same* variant set and any difference between them is method, not data.
#   * No absolute paths. Callers build paths with pp() (see below) plus the config.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

color_scheme <- c("#00274c", "#ffcb05")   # maize & blue, matching the rest of the project

# ---- project paths ------------------------------------------------------------
# Deliberately NOT here::here(). This directory has no .here/.Rproj/.git marker, so
# here() anchors on the repo root one level up and every relative path misses. We
# locate the root by walking up for a sentinel file instead, which works whether
# the working directory is the document folder, the analysis folder, or a cluster
# scratch copy — and does not disturb path resolution for the alzheimers/ scripts.

find_project_root <- function(anchor = "config_robust_mr.yml", start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, anchor))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop("Could not find '", anchor, "' in ", start, " or any parent directory.\n",
           "  Set it explicitly:  set_project_root('/path/to/Human Genetics')")
    }
    d <- parent
  }
}

set_project_root <- function(path = find_project_root()) {
  options(robust_mr.root = normalizePath(path, mustWork = TRUE))
  invisible(getOption("robust_mr.root"))
}

#' Build a path relative to the project root. Drop-in replacement for here::here().
#' Absolute paths (and ~) pass through untouched, so config entries can point at
#' shared cluster locations — e.g. an LD panel already on /nfs/turbo — instead of
#' being forced under the project root.
pp <- function(...) {
  parts <- as.character(c(...))
  if (length(parts) && grepl("^(/|~)", parts[1])) {
    return(path.expand(do.call(file.path, as.list(parts))))
  }
  root <- getOption("robust_mr.root")
  if (is.null(root)) root <- set_project_root()
  file.path(root, ...)
}

# ---- config -------------------------------------------------------------------

load_robust_cfg <- function(path = pp("config_robust_mr.yml")) {
  yaml::read_yaml(path)
}

#' Resolve the plink binary, in precedence order:
#'   1. explicit `plink_bin` argument
#'   2. the PLINK_BIN environment variable  <- how the Great Lakes sbatch passes
#'      the module-loaded plink, since genetics.binaRies may not be installed there
#'   3. paths$plink_bin in config_robust_mr.yml
#'   4. the genetics.binaRies bundled copy (what ieugwasr/TwoSampleMR use locally)
#'   5. plink on $PATH
resolve_plink <- function(plink_bin = NULL, cfg = NULL) {
  ok <- function(x) !is.null(x) && length(x) == 1 && nzchar(x) && !identical(x, "NULL")

  if (ok(plink_bin)) return(plink_bin)

  env <- Sys.getenv("PLINK_BIN", unset = "")
  if (ok(env)) return(env)

  if (!is.null(cfg) && ok(cfg$paths$plink_bin)) return(cfg$paths$plink_bin)

  if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
    b <- try(genetics.binaRies::get_plink_binary(), silent = TRUE)
    if (!inherits(b, "try-error") && ok(b)) return(b)
  }

  for (nm in c("plink", "plink1.9", "plink19")) {
    onpath <- unname(Sys.which(nm))
    if (ok(onpath)) return(onpath)
  }

  # plink2 is NOT a substitute: ieugwasr::ld_clump() calls --clump with plink 1.9
  # syntax and output conventions.
  hint <- if (nzchar(Sys.which("plink2")))
    "\n  NOTE: plink2 is on your PATH but will not work — ld_clump() needs plink 1.9." else ""

  stop("No plink 1.9 binary found. Easiest fix:\n",
       "    remotes::install_github('MRCIEU/genetics.binaRies')\n",
       "  Or set the PLINK_BIN environment variable, or paths$plink_bin in ",
       "config_robust_mr.yml.", hint)
}

# ---- dependency management ----------------------------------------------------

#' Check (and optionally install) the packages this arm needs.
#' CAUSE and MR-APSS are GitHub-only; both pull a non-trivial dependency tree.
check_robust_mr_deps <- function(install = FALSE) {
  cran   <- c("tidyverse", "data.table", "yaml", "here", "R.utils",
              "ieugwasr", "TwoSampleMR", "knitr", "kableExtra")
  github <- c(cause = "jean997/cause", MRAPSS = "YangLabHKUST/MR-APSS",
              mixsqp = "stephenslab/mixsqp", ashr = "stephens999/ashr")

  status <- tibble::tibble(
    package = c(cran, names(github)),
    source  = c(rep("CRAN", length(cran)), unname(github)),
    installed = purrr::map_lgl(c(cran, names(github)),
                               ~ requireNamespace(.x, quietly = TRUE))
  )

  if (install && any(!status$installed)) {
    miss <- status %>% dplyr::filter(!installed)
    for (i in seq_len(nrow(miss))) {
      if (miss$source[i] == "CRAN") {
        install.packages(miss$package[i])
      } else {
        remotes::install_github(miss$source[i])
      }
    }
    status$installed <- purrr::map_lgl(status$package,
                                       ~ requireNamespace(.x, quietly = TRUE))
  }
  status
}

# ---- reading ------------------------------------------------------------------

#' Print the header of a (possibly gzipped) summary statistics file.
#' Use this before filling in the `cols:` block of config_robust_mr.yml — the exact
#' column names in GLGC / GEFOS releases have changed between versions and guessing
#' them is how silent allele flips happen.
inspect_gwas_header <- function(path, n = 3) {
  stopifnot(file.exists(path))
  con <- if (grepl("\\.gz$|\\.bgz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  readLines(con, n = n)
}

#' Read a GWAS summary statistics file into canonical columns.
#'
#' @param path  file path (plain or gzipped)
#' @param cols  named list mapping canonical -> file column name. Recognised
#'   canonical names: snp, chr, bp, a1 (EFFECT allele), a2, beta, se, pval, eaf,
#'   n, z, info. Supply either (beta, se) or z; either n or `n_fixed`.
#' @param n_fixed  scalar N to use when the file has no per-SNP N column.
#' @param log10p  TRUE if the p-value column is -log10(p) (GLGC 2021 ships both).
read_gwas <- function(path, cols, n_fixed = NULL, log10p = FALSE) {
  stopifnot(file.exists(path))

  # Only request columns that actually exist. fread() errors on a missing `select`
  # name, and upstream headers change between releases, so intersect first and say
  # loudly what was dropped rather than failing on the whole file.
  available <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  ren  <- setNames(unlist(cols, use.names = FALSE), names(cols))
  gone <- ren[!ren %in% available]
  if (length(gone)) {
    warning("Columns in the config that are not in ", basename(path), ": ",
            paste(sprintf("%s -> '%s'", names(gone), gone), collapse = ", "),
            "\n  Header is: ", paste(available, collapse = ", "),
            call. = FALSE)
    ren <- ren[ren %in% available]
  }
  if (!all(c("snp", "a1", "a2") %in% names(ren))) {
    stop("Need at least snp, a1 and a2 columns mapped for ", basename(path))
  }

  dat <- data.table::fread(path, select = unname(ren), showProgress = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::rename(dplyr::all_of(ren))

  dat <- dat %>%
    dplyr::mutate(
      snp = as.character(.data$snp),
      a1  = toupper(as.character(.data$a1)),
      a2  = toupper(as.character(.data$a2))
    )

  if (!"n" %in% names(dat)) {
    if (is.null(n_fixed)) stop("No N column and no n_fixed supplied for ", basename(path))
    dat$n <- n_fixed
  }
  if (isTRUE(log10p) && "pval" %in% names(dat))  dat$pval <- 10^(-dat$pval)

  # derive z from beta/se, or beta/se from z + n, whichever is missing
  if (!"z" %in% names(dat) && all(c("beta", "se") %in% names(dat))) {
    dat$z <- dat$beta / dat$se
  }
  if (!"z" %in% names(dat) && !"beta" %in% names(dat)) {
    stop("Need either (beta, se) or z mapped for ", basename(path))
  }
  if (!"beta" %in% names(dat) && "z" %in% names(dat)) {
    # standardized scale: b = Z/sqrt(N), se = 1/sqrt(N) (the MR-APSS convention)
    dat$beta <- dat$z / sqrt(dat$n)
    dat$se   <- 1 / sqrt(dat$n)
  }
  if (!"pval" %in% names(dat)) {
    dat$pval <- 2 * pnorm(abs(dat$z), lower.tail = FALSE)
  } else {
    check_pval_consistency(dat, basename(path))
  }
  dat
}

#' Warn if a file's p-value column disagrees with the one implied by beta/se.
#'
#' Catches mislabelled or shifted columns, which are silent otherwise. Real case:
#' the GWAS Catalog reformatting of Morris 2019 (GCST006979) labels one BOLT-LMM
#' p-value column `p_value` and the other `n`, so the column named `p_value` does
#' not correspond to the reported beta/se.
check_pval_consistency <- function(dat, label, n_check = 20000, tol = 0.05) {
  ok <- which(is.finite(dat$pval) & dat$pval > 1e-300 & dat$pval < 1 &
              is.finite(dat$z) & abs(dat$z) > 0.1)
  if (length(ok) < 100) return(invisible(NULL))
  idx <- if (length(ok) > n_check) sample(ok, n_check) else ok

  derived <- 2 * pnorm(abs(dat$z[idx]), lower.tail = FALSE)
  lr <- stats::median(log10(dat$pval[idx] / derived), na.rm = TRUE)

  if (!is.finite(lr) || abs(lr) > tol) {
    warning("[", label, "] the p-value column disagrees with beta/se ",
            "(median log10 ratio = ", round(lr, 3), ").\n",
            "  The mapped p-value column is probably not the one matching these ",
            "effect estimates.\n",
            "  Safest fix: remove `pval:` from this file's cols: block and let ",
            "read_gwas() derive it from beta/se.", call. = FALSE)
  } else {
    message("[", label, "] p-value column is consistent with beta/se ",
            "(median log10 ratio = ", round(lr, 3), ").")
  }
  invisible(lr)
}

#' Read an OpenGWAS GWAS-VCF file into the same canonical columns as read_gwas().
#'
#' Use this instead of read_gwas() when you would rather pull the outcome from
#' OpenGWAS (https://gwas.mrcieu.ac.uk/files/<id>/<id>.vcf.gz) than chase a
#' publisher's flat file. NOTE: this is the *bulk file*, not the API —
#' ieugwasr::associations() / extract_outcome_data() are per-SNP queries and are
#' not a sensible way to retrieve a genome-wide dataset.
#'
#' GWAS-VCF puts the stats in FORMAT fields: ES (effect size), SE, LP (-log10 p),
#' AF (effect-allele freq), SS (sample size). ALT is the effect allele.
#'
#' Prefers bcftools (fast, low memory). Falls back to VariantAnnotation, which
#' reads the whole file into memory and will struggle above a few million rows.
read_gwasvcf <- function(path, n_fixed = NULL) {
  stopifnot(file.exists(path))
  bcftools <- unname(Sys.which("bcftools"))

  if (nzchar(bcftools)) {
    fmt <- "%ID\t%CHROM\t%POS\t%ALT\t%REF\t[%ES]\t[%SE]\t[%LP]\t[%AF]\t[%SS]\n"
    tsv <- tempfile(fileext = ".tsv")
    on.exit(unlink(tsv), add = TRUE)
    st <- system2(bcftools, c("query", "-f", shQuote(fmt), shQuote(path)), stdout = tsv)
    if (st != 0) stop("bcftools query failed on ", basename(path))
    dat <- data.table::fread(tsv, na.strings = c(".", "NA", ""),
                             col.names = c("snp", "chr", "bp", "a1", "a2",
                                           "beta", "se", "lp", "eaf", "n"),
                             showProgress = FALSE) %>%
      tibble::as_tibble()
  } else {
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
      stop("Reading GWAS-VCF needs either bcftools on $PATH (brew install bcftools) ",
           "or the VariantAnnotation Bioconductor package.")
    }
    message("bcftools not found; falling back to VariantAnnotation (memory heavy).")
    v   <- VariantAnnotation::readVcf(path)
    geno <- VariantAnnotation::geno(v)
    g <- function(f) if (!is.null(geno[[f]])) as.numeric(geno[[f]][, 1]) else NA_real_
    rr <- SummarizedExperiment::rowRanges(v)
    dat <- tibble::tibble(
      snp  = names(rr),
      chr  = as.character(GenomicRanges::seqnames(rr)),
      bp   = GenomicRanges::start(rr),
      a1   = as.character(unlist(VariantAnnotation::alt(v))),
      a2   = as.character(VariantAnnotation::ref(v)),
      beta = g("ES"), se = g("SE"), lp = g("LP"), eaf = g("AF"), n = g("SS")
    )
  }

  dat <- dat %>%
    dplyr::mutate(
      snp  = as.character(.data$snp),
      chr  = suppressWarnings(as.integer(gsub("^chr", "", .data$chr))),
      a1   = toupper(as.character(.data$a1)),
      a2   = toupper(as.character(.data$a2)),
      pval = 10^(-.data$lp),                 # GWAS-VCF stores -log10(p) as LP
      z    = .data$beta / .data$se
    ) %>%
    dplyr::select(-"lp")

  if (all(is.na(dat$n))) {
    if (is.null(n_fixed)) stop("No SS field in ", basename(path), " and no n_fixed given")
    dat$n <- n_fixed
  }
  dat %>% dplyr::filter(!is.na(.data$beta), !is.na(.data$se), .data$se > 0)
}

# ---- QC -----------------------------------------------------------------------

AMBIGUOUS <- c("AT", "TA", "CG", "GC")
VALID_ALLELES <- c("A", "C", "G", "T")

#' MR-APSS-style QC. Applied identically to both traits so CAUSE and MR-APSS
#' operate on the same variant set.
#'
#' @param dat canonical GWAS tibble from read_gwas()
#' @param hm3 character vector of HapMap3 rsIDs (from w_hm3.snplist), or NULL to skip
#' @param maf_min,info_min thresholds; skipped silently if the column is absent
#' @param drop_mhc exclude chr6:26-34Mb (GRCh37)
qc_gwas <- function(dat, hm3 = NULL, maf_min = 0.05, info_min = 0.9,
                    drop_mhc = TRUE, chi2_max = NULL, label = "trait") {

  n0 <- nrow(dat)
  log <- tibble::tibble(step = "input", n = n0)
  add <- function(log, step, d) dplyr::bind_rows(log, tibble::tibble(step = step, n = nrow(d)))

  dat <- dat %>% dplyr::filter(!is.na(.data$snp), !duplicated(.data$snp))
  log <- add(log, "unique rsID", dat)

  dat <- dat %>% dplyr::filter(.data$a1 %in% VALID_ALLELES, .data$a2 %in% VALID_ALLELES)
  log <- add(log, "biallelic ACGT", dat)

  dat <- dat %>% dplyr::filter(!paste0(.data$a1, .data$a2) %in% AMBIGUOUS, .data$a1 != .data$a2)
  log <- add(log, "unambiguous strand", dat)

  if (!is.null(hm3)) {
    dat <- dat %>% dplyr::filter(.data$snp %in% hm3)
    log <- add(log, "in HapMap3", dat)
  }
  if ("eaf" %in% names(dat) && !is.na(maf_min)) {
    dat <- dat %>% dplyr::filter(pmin(.data$eaf, 1 - .data$eaf) > maf_min)
    log <- add(log, paste0("MAF > ", maf_min), dat)
  }
  if ("info" %in% names(dat) && !is.na(info_min)) {
    dat <- dat %>% dplyr::filter(.data$info > info_min)
    log <- add(log, paste0("INFO > ", info_min), dat)
  }
  if (drop_mhc) {
    if (all(c("chr", "bp") %in% names(dat))) {
      dat <- dat %>% dplyr::filter(!(.data$chr == 6 & .data$bp > 26e6 & .data$bp < 34e6))
      log <- add(log, "MHC removed", dat)
    } else {
      message("[", label, "] no chr/bp mapped, so MHC could not be excluded here. ",
              "This is usually harmless: the merge with a dataset that DOES have ",
              "coordinates (the outcome) removes them from the analysis anyway. ",
              "Map chr/bp in config_robust_mr.yml if you want it done symmetrically.")
    }
  }

  # chi2 cap: removes a handful of implausibly large statistics that destabilise LDSC
  if (is.null(chi2_max)) chi2_max <- max(stats::median(dat$n, na.rm = TRUE) / 1000, 80)
  dat <- dat %>% dplyr::filter(.data$z^2 < chi2_max)
  log <- add(log, sprintf("chi2 < %.0f", chi2_max), dat)

  attr(dat, "qc_log")   <- log %>% dplyr::mutate(trait = label, .before = 1)
  attr(dat, "chi2_max") <- chi2_max
  dat
}

qc_log <- function(dat) attr(dat, "qc_log")

#' Read per-chromosome LDSC LD scores, tolerating either distribution layout.
#'
#' MRAPSS::est_paras(ldscore.dir=) hard-codes the `eur_w_ld_chr` naming
#' (`1.l2.ldscore.gz`, `1.l2.M_5_50`). The S-LDSC bundles on Zenodo instead ship
#' `LDscore.1.l2.ldscore.gz`, so pointing est_paras() at them fails. This reads
#' either layout and returns the `ld` data.frame and `M` scalar that est_paras()
#' also accepts directly — so it does not matter which bundle you managed to get.
#'
#' @return list(ld = data.frame with SNP and L2, M = total SNPs in the LD score calc)
read_ldscores <- function(dir) {
  if (!dir.exists(dir)) stop("LD score directory not found: ", dir)

  # recursive: the S-LDSC tarballs unpack into a nested LDscore/ subdirectory, and
  # the anchored "(^|\\.)" prefix keeps chr 1 from also matching chr 11/21.
  pick <- function(chr, suffix) {
    for (pat in c(sprintf("(^|\\.)%d\\.l2\\.%s$", chr, suffix),
                  sprintf("(^|\\.)%d\\.l2\\.%s\\.gz$", chr, suffix))) {
      hits <- list.files(dir, pattern = pat, full.names = TRUE, recursive = TRUE)
      if (length(hits)) return(hits[1])
    }
    NA_character_
  }

  ld <- purrr::map_dfr(1:22, function(chr) {
    f <- pick(chr, "ldscore")
    if (is.na(f)) {
      found <- list.files(dir, recursive = TRUE)
      stop("No LD score file for chr", chr, " under ", dir,
           "\n  Expected 1.l2.ldscore.gz or LDscore.1.l2.ldscore.gz",
           "\n  Did you untar the archive? Found ", length(found), " file(s): ",
           paste(utils::head(found, 5), collapse = ", "))
    }
    data.table::fread(f, showProgress = FALSE) %>% tibble::as_tibble()
  })

  m_files <- purrr::map_chr(1:22, ~ pick(.x, "M_5_50"))
  M <- if (any(is.na(m_files))) {
    warning("No .M_5_50 files found; falling back to M = nrow(ld). ",
            "This changes the scale of Omega but not the causal estimate.")
    nrow(ld)
  } else {
    sum(purrr::map_dbl(m_files, ~ sum(as.numeric(readLines(.x)))))
  }

  if (!"L2" %in% names(ld)) {
    stop("LD score files have no L2 column; got: ", paste(names(ld), collapse = ", "))
  }
  message(sprintf("Read LD scores for %s SNPs (M = %s)",
                  format(nrow(ld), big.mark = ","), format(M, big.mark = ",")))
  list(ld = as.data.frame(ld), M = M)
}

#' Read the LDSC HapMap3 snplist (w_hm3.snplist).
read_hm3 <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    warning("HapMap3 snplist not found; skipping the HM3 restriction. ",
            "LDSC-based estimates (MR-APSS C and Omega) are not reliable without it.")
    return(NULL)
  }
  data.table::fread(path, showProgress = FALSE)$SNP
}

# ---- format converters --------------------------------------------------------

#' Canonical -> MR-APSS input: SNP, A1, A2, Z, P, N, chi2.
#'
#' chi2 is NOT optional even though it is trivially Z^2. est_paras() merges dat1
#' and dat2 by SNP and hands the result to ldsc_GC(), which selects
#' c("SNP","chi2.x","chi2.y","N.x","N.y","Zxy","L2") — it derives Zxy itself but
#' expects chi2 to already be present on each side. Omitting it fails deep in the
#' call stack with "undefined columns selected". MR-APSS's own format_data()
#' adds it; we add it here so we can keep control of the QC.
to_apss <- function(dat) {
  dat %>%
    dplyr::transmute(SNP = .data$snp, A1 = .data$a1, A2 = .data$a2,
                     Z = .data$z, P = .data$pval, N = .data$n,
                     chi2 = .data$z^2) %>%
    dplyr::filter(is.finite(.data$Z), is.finite(.data$N))
}

#' Canonical -> CAUSE input for gwas_merge().
to_cause <- function(dat, standardize = TRUE) {
  out <- dat %>%
    dplyr::transmute(SNP = .data$snp, A1 = .data$a1, A2 = .data$a2,
                     beta_hat = .data$beta, se = .data$se, pval = .data$pval,
                     N = .data$n)
  if (standardize) {
    # Put both traits on a per-SD scale so gamma is directly comparable to the
    # MR-APSS beta (which is always Z/sqrt(N) internally). Harmless if the source
    # betas are already per-SD, essential if they are not.
    out <- out %>% dplyr::mutate(beta_hat = .data$beta_hat / (.data$se * sqrt(.data$N)),
                                 se = 1 / sqrt(.data$N))
  }
  out %>% dplyr::filter(is.finite(.data$beta_hat), is.finite(.data$se), .data$se > 0)
}

# ---- LD pruning / clumping ----------------------------------------------------

#' LD-prune a set of variants with local plink via ieugwasr.
#' Used for both CAUSE (r2 0.01-0.1, p<1e-3) and MR-APSS (r2 0.001, p<5e-5).
clump_local <- function(dat, snp_col = "SNP", p_col = "pval",
                        r2 = 0.01, kb = 10000, p_thresh = 1e-3,
                        bfile = NULL, plink_bin = NULL) {
  if (!is.null(bfile)) plink_bin <- resolve_plink(plink_bin)
  keep <- dat %>% dplyr::filter(.data[[p_col]] < p_thresh)
  if (nrow(keep) == 0) return(keep[0, ])
  message(sprintf("Clumping %s variants at p<%.0e, r2=%.3f, kb=%d ...",
                  format(nrow(keep), big.mark = ","), p_thresh, r2, kb))
  # as.character(): est_paras()/clump() can hand back factor SNP columns, which
  # ld_clump writes to disk as integer codes.
  cl <- ieugwasr::ld_clump(
    dplyr::tibble(rsid = as.character(keep[[snp_col]]),
                  pval = as.numeric(keep[[p_col]])),
    clump_r2 = r2, clump_kb = kb, clump_p = p_thresh,
    bfile = bfile, plink_bin = plink_bin
  )
  keep %>% dplyr::filter(as.character(.data[[snp_col]]) %in% cl$rsid)
}

#' Make a clumped MR-APSS dataset safe to pass to MRAPSS().
#'
#' MRAPSS() reads its selection threshold as `unique(MRdat$Threshold)` and feeds it
#' straight to qnorm(Threshold/2) with no coercion, so a Threshold column that is
#' character/factor (as some clump() versions produce) fails deep in the EM with
#' "non-numeric argument to binary operator". The threshold is ours to define, so
#' set it explicitly rather than inheriting whatever clump() wrote. Also coerces
#' the numeric columns the EM relies on, and reports anything it had to change.
sanitise_mrdat <- function(MRdat, iv_threshold) {
  MRdat <- as.data.frame(MRdat)
  notes <- character()

  before <- if ("Threshold" %in% names(MRdat)) class(MRdat$Threshold)[1] else "absent"
  if (!identical(before, "numeric")) {
    notes <- c(notes, sprintf("Threshold was %s -> set to %.0e", before, iv_threshold))
  }
  MRdat$Threshold <- as.numeric(iv_threshold)

  num_cols <- c("b.exp", "b.out", "se.exp", "se.out", "pval.exp", "pval.out", "L2")
  for (nm in intersect(num_cols, names(MRdat))) {
    if (!is.numeric(MRdat[[nm]])) {
      notes <- c(notes, sprintf("%s was %s -> as.numeric", nm, class(MRdat[[nm]])[1]))
      MRdat[[nm]] <- as.numeric(as.character(MRdat[[nm]]))
    }
  }

  need <- intersect(c("b.exp", "b.out", "se.exp", "se.out", "L2"), names(MRdat))
  keep <- stats::complete.cases(MRdat[, need, drop = FALSE])
  if (any(!keep)) notes <- c(notes, sprintf("dropped %d row(s) with NA", sum(!keep)))
  MRdat <- MRdat[keep, , drop = FALSE]

  if (length(notes)) message("  sanitise_mrdat: ", paste(notes, collapse = "; "))
  if (nrow(MRdat) < 4) stop("Fewer than 4 instruments after sanitising — MRAPSS needs >= 4.")
  MRdat
}

# ---- result tidiers -----------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Pull a numeric scalar out of a result list, coercing from character if needed.
#' MRAPSS returns `pvalue` as a CHARACTER string ("1.8746e-16"), so a plain
#' is.numeric() guard silently yields NA. Verified against MRAPSS on the bundled
#' BMI->T2D example, 2026-08-10.
pluck_num <- function(x, ...) {
  for (nm in c(...)) {
    v <- x[[nm]]
    if (!is.null(v) && length(v) >= 1) {
      out <- suppressWarnings(as.numeric(v[1]))
      if (!is.na(out)) return(out)
    }
  }
  NA_real_
}

#' MR-APSS result -> one-row tibble.
#'
#' Field names verified against the installed package via str() on the bundled
#' example: the object carries beta, beta.se, pvalue (character), tau.sq, sigma.sq,
#' pi0, IVsignal.sum, Threshold, method, MRdat, and post (mu, Pi, IVsignal.sum).
#' There is no instrument-count field — it is nrow(MRdat).
tidy_apss <- function(res, exposure, outcome, arm, C_setting = "corrected") {
  b  <- pluck_num(res, "beta")
  se <- pluck_num(res, "beta.se")

  tibble::tibble(
    method    = "MR-APSS",
    arm       = arm,
    C_setting = C_setting,
    exposure  = exposure,
    outcome   = outcome,
    b = b, se = se,
    lci = b - 1.96 * se, uci = b + 1.96 * se,
    pval = pluck_num(res, "pvalue"),
    n_iv = if (!is.null(res$MRdat)) nrow(res$MRdat) else NA_integer_,
    # Posterior-expected count of instruments carrying foreground (valid) signal:
    # post$Pi is the per-SNP posterior P(Z_j = 1), so its sum is the effective
    # number of valid IVs. Reported alongside pi0 rather than derived from it,
    # because the pi0 -> "valid IV" mapping is not documented unambiguously.
    n_valid_iv = if (!is.null(res$post$Pi)) sum(as.numeric(res$post$Pi)) else NA_real_,
    pi0        = pluck_num(res, "pi0"),
    iv_signal  = pluck_num(res, "IVsignal.sum"),
    threshold  = pluck_num(res, "Threshold"),
    tau2       = pluck_num(res, "tau.sq"),
    sigma2     = pluck_num(res, "sigma.sq")
  )
}

#' CAUSE result -> one-row tibble. `gamma` is the causal effect under the causal
#' model; the *test* is the sharing-vs-causal ELPD comparison, not a Wald p on gamma.
tidy_cause <- function(res, exposure, outcome, arm, ci_size = 0.95) {
  s  <- summary(res, ci_size = ci_size)
  # quants[[2]] = causal model: columns gamma, eta, q; rows median, lower, upper
  q  <- s$quants[[2]]
  el <- as.data.frame(res$elpd)
  row <- el %>% dplyr::filter(.data$model1 == "sharing", .data$model2 == "causal")
  if (nrow(row) == 0) row <- el[nrow(el), ]

  tibble::tibble(
    method   = "CAUSE",
    arm      = arm,
    exposure = exposure,
    outcome  = outcome,
    b   = as.numeric(q[1, "gamma"]),
    lci = as.numeric(q[2, "gamma"]),
    uci = as.numeric(q[3, "gamma"]),
    eta_med = as.numeric(q[1, "eta"]),
    q_med   = as.numeric(q[1, "q"]),
    delta_elpd    = as.numeric(row$delta_elpd),
    se_delta_elpd = as.numeric(row$se_delta_elpd),
    z    = as.numeric(row$z),
    pval = pnorm(as.numeric(row$z), lower.tail = TRUE),
    verdict = dplyr::case_when(
      pnorm(as.numeric(row$z), lower.tail = TRUE) < 0.05 ~ "causal preferred over sharing",
      TRUE ~ "sharing not rejected"
    )
  )
}

# ---- Bayesian reading of CAUSE ------------------------------------------------
# CAUSE fits full posteriors for gamma (causal effect), eta (shared-factor effect)
# and q (proportion of variants acting through the shared factor). The ELPD z-test
# is a frequentist wrapper on top of that. These helpers go back to the posterior.

#' Posterior quantiles for a CAUSE model's parameters.
#'
#' summary.cause(ci_size = c) returns the median plus the (1-c)/2 and (1+c)/2
#' quantiles. Sweeping ci_size therefore recovers an arbitrary quantile ladder
#' without touching CAUSE's internals — which is the robust way to do this, since
#' the layout of the fitted object has changed between versions.
cause_posterior_quantiles <- function(fit, model = c("causal", "sharing"),
                                      probs = c(0.025, 0.05, 0.1, 0.25,
                                                0.75, 0.9, 0.95, 0.975)) {
  model <- match.arg(model)
  idx   <- if (model == "causal") 2L else 1L

  med    <- summary(fit, ci_size = 0.5)$quants[[idx]]
  params <- colnames(med)

  rows <- purrr::map_dfr(sort(unique(probs)), function(p) {
    s   <- summary(fit, ci_size = abs(1 - 2 * p))$quants[[idx]]
    row <- if (p < 0.5) 2L else 3L
    tibble::tibble(param = params, prob = p, value = as.numeric(s[row, ]))
  })
  rows <- dplyr::bind_rows(
    rows,
    tibble::tibble(param = params, prob = 0.5, value = as.numeric(med[1, ]))
  )

  rows %>%
    dplyr::arrange(.data$param, .data$prob) %>%
    tidyr::pivot_wider(names_from = "prob", values_from = "value",
                       names_prefix = "q") %>%
    dplyr::mutate(model = model, .before = 1)
}

#' Posterior probability that a parameter lies below a threshold, interpolated
#' from the quantile ladder. P(gamma < 0) is the natural "is the effect negative"
#' statement, and needs no null hypothesis.
cause_p_below <- function(qtab, param = "gamma", threshold = 0) {
  row <- qtab %>% dplyr::filter(.data$param == !!param)
  if (nrow(row) == 0) return(NA_real_)
  qcols <- grep("^q0", names(row), value = TRUE)
  qs <- as.numeric(unlist(row[, qcols]))
  ps <- as.numeric(sub("^q", "", qcols))
  o  <- order(qs)
  stats::approx(qs[o], ps[o], xout = threshold, rule = 2)$y
}

#' Pseudo-BMA model weights from the CAUSE ELPD table.
#'
#' CAUSE reports delta_elpd = elpd(model1) - elpd(model2). Anchoring on the null
#' model, weight_k proportional to exp(elpd_k) gives an interpretable share of
#' predictive support across null / sharing / causal — a far better summary of an
#' inconclusive comparison than a p-value on one pairwise contrast.
#'
#' NOTE: these are naive (pseudo-BMA) weights. They ignore the standard error on
#' delta_elpd, which here is comparable to delta_elpd itself, so treat them as a
#' point summary of a very uncertain quantity. cause_model_weights_boot() gives
#' the uncertainty.
cause_model_weights <- function(elpd_tab) {
  g <- function(m1, m2) {
    v <- elpd_tab$delta_elpd[elpd_tab$model1 == m1 & elpd_tab$model2 == m2]
    if (length(v) == 0) NA_real_ else v[1]
  }
  elpd <- c(null = 0, sharing = -g("null", "sharing"), causal = -g("null", "causal"))
  w <- exp(elpd - max(elpd)); w <- w / sum(w)
  tibble::tibble(model = names(elpd), elpd_rel_null = as.numeric(elpd),
                 weight = as.numeric(w))
}

#' Uncertainty on the model weights, by resampling the ELPD differences from
#' their reported standard errors. Shows how much of the "74% causal" is signal.
cause_model_weights_boot <- function(elpd_tab, n = 10000, seed = 1) {
  set.seed(seed)
  g <- function(m1, m2, col) {
    v <- elpd_tab[[col]][elpd_tab$model1 == m1 & elpd_tab$model2 == m2]
    if (length(v) == 0) NA_real_ else v[1]
  }
  draws <- purrr::map_dfr(seq_len(n), function(i) {
    e_sh <- rnorm(1, -g("null", "sharing", "delta_elpd"),
                  g("null", "sharing", "se_delta_elpd"))
    e_ca <- rnorm(1, -g("null", "causal", "delta_elpd"),
                  g("null", "causal", "se_delta_elpd"))
    e <- c(null = 0, sharing = e_sh, causal = e_ca)
    w <- exp(e - max(e)); w <- w / sum(w)
    tibble::tibble(model = names(w), weight = as.numeric(w))
  })
  draws %>%
    dplyr::group_by(.data$model) %>%
    dplyr::summarise(weight_med = stats::median(.data$weight),
                     weight_lo  = stats::quantile(.data$weight, 0.025),
                     weight_hi  = stats::quantile(.data$weight, 0.975),
                     .groups = "drop")
}

#' Summarise the MR-APSS background parameters. C[1,2] is the cross-trait LDSC
#' intercept: a *direct empirical estimate of sample overlap* (plus any shared
#' stratification). It should be ~0 for genuinely independent samples.
tidy_background <- function(paras, exposure, outcome, arm) {
  C <- paras$C; O <- paras$Omega
  tibble::tibble(
    arm = arm, exposure = exposure, outcome = outcome,
    C11 = C[1, 1], C22 = C[2, 2], C12 = C[1, 2],
    overlap_flag = abs(C[1, 2]) > 0.02,
    Omega11 = O[1, 1], Omega22 = O[2, 2], Omega12 = O[1, 2],
    rg = O[1, 2] / sqrt(O[1, 1] * O[2, 2])
  )
}

# ---- plotting -----------------------------------------------------------------

#' Forest plot comparing robust methods against the conventional estimates.
robust_mr_forest <- function(dat, title = NULL,
                             xlab = "Effect on eBMD (SD per SD higher LDL-C)") {
  dat %>%
    dplyr::mutate(label = paste0(.data$method,
                                 dplyr::if_else(is.na(.data$arm), "", paste0(" — ", .data$arm))),
                  label = forcats::fct_rev(forcats::fct_inorder(.data$label))) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$b, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lci, xmax = .data$uci),
                            height = 0.18, colour = color_scheme[1]) +
    ggplot2::geom_point(size = 2.6, colour = color_scheme[1]) +
    ggplot2::labs(x = xlab, y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 13)
}
