# helpers.R — shared utilities for the smoking->AD nAChR decomposition pipeline.
# Sourced by each layer's .qmd. Keep functions pure and side-effect free except
# where noted (writers). No absolute paths: callers pass paths via here::here().

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# Maize & blue, matching the rest of the project.
color_scheme <- c("#00274c", "#ffcb05")

# Standard error of the mean.
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x))

# Read the project config once.
load_config <- function() yaml::read_yaml(here::here("config.yml"))

# ---- Gene annotation ---------------------------------------------------------
# Annotate a data frame of SNPs (cols: SNP, chr, pos in hg19/GRCh37) with the
# nearest gene symbol + distance. Returns input joined with nearest_gene,distance.
annotate_nearest_gene <- function(df, chr_col = "chr.exposure", pos_col = "pos.exposure",
                                  snp_col = "SNP") {
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
  })
  snps <- GRanges(
    seqnames = paste0("chr", df[[chr_col]]),
    ranges   = IRanges(start = as.numeric(df[[pos_col]]),
                       end   = as.numeric(df[[pos_col]])),
    rsid     = df[[snp_col]]
  )
  txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- suppressWarnings(genes(txdb))
  nidx  <- nearest(snps, genes)
  dist  <- distance(snps, genes[nidx])
  syms  <- mapIds(org.Hs.eg.db, keys = genes[nidx]$gene_id,
                  column = "SYMBOL", keytype = "ENTREZID")
  ann <- tibble(!!snp_col := snps$rsid,
                nearest_gene = unname(syms),
                distance = dist)
  dplyr::left_join(df, ann, by = snp_col)
}

# ---- Mechanism binning -------------------------------------------------------
# GRCh37/hg19 coordinates of the panel genes, with their mechanism bin. Used to bin SNPs by
# cis-window membership rather than literal nearest gene — so e.g. the 15q25 cluster lead
# (technically nearest HYKK) is correctly assigned to the nAChR bin, and the chr19 CYP2A6
# locus (technically nearest a pseudogene/EGLN2) to nicotine_metabolism.
panel_gene_coords <- function() {
  tibble::tribble(
    ~gene,     ~bin,                    ~chr,       ~start,      ~end,
    "CHRNA5",  "nAChR_pharmacodynamic",   15L,  78857862,  78885393,
    "CHRNA3",  "nAChR_pharmacodynamic",   15L,  78885394,  78919647,
    "CHRNB4",  "nAChR_pharmacodynamic",   15L,  78919290,  78937340,
    "CHRNA4",  "nAChR_pharmacodynamic",   20L,  61975397,  62019864,
    "CHRNB2",  "nAChR_pharmacodynamic",    1L, 154543992, 154579816,
    "CHRNA7",  "nAChR_pharmacodynamic",   15L,  32322677,  32464722,
    "CHRNA6",  "nAChR_pharmacodynamic",    8L,  42634258,  42648862,
    "CHRNB3",  "nAChR_pharmacodynamic",    8L,  42551452,  42591443,
    "CYP2A6",  "nicotine_metabolism",     19L,  41349443,  41356352,
    "CYP2B6",  "nicotine_metabolism",     19L,  41497203,  41524301,
    "DRD2",    "reward_behavioral",       11L, 113280318, 113346413,
    "DBH",     "reward_behavioral",        9L, 136501486, 136524466)
}

# Assign each SNP a mechanism bin by cis-window membership of the panel genes (flank =
# cfg$cis_window_kb). A SNP within a panel gene's window is assigned that gene's bin (nearest
# panel gene wins ties); otherwise "other_pleiotropic". Adds `bin_gene` (the panel gene that
# drove the call) for transparency. Falls back to nearest-gene symbol matching only if the
# data frame lacks position columns.
assign_mechanism_bin <- function(df, cfg, chr_col = "chr.exposure", pos_col = "pos.exposure",
                                 gene_col = "nearest_gene") {
  coords <- panel_gene_coords()
  flank <- cfg$cis_window_kb * 1000

  if (!all(c(chr_col, pos_col) %in% names(df))) {
    # Fallback: exact nearest-gene symbol membership (legacy behaviour).
    panel <- cfg$gene_panel
    bin_of <- function(g) {
      if (is.na(g)) return("other_pleiotropic")
      for (b in names(panel)) if (g %in% panel[[b]]) return(b)
      "other_pleiotropic"
    }
    return(df %>% mutate(mechanism = vapply(.data[[gene_col]], bin_of, character(1)),
                         bin_gene = NA_character_))
  }

  # Hybrid rule: (1) cis-window membership of a panel gene, else (2) nearest gene is itself a
  # panel gene (catches loci just outside the window, e.g. DRD2). Window takes precedence.
  panel_sym <- coords %>% dplyr::select(gene, bin)
  assign_one <- function(chr, pos, ng) {
    chr <- suppressWarnings(as.integer(chr)); pos <- as.numeric(pos)
    hits <- coords %>%
      dplyr::filter(chr == !!chr, pos >= start - flank, pos <= end + flank)
    if (nrow(hits) > 0) {
      hits <- hits %>%
        dplyr::mutate(dist = ifelse(pos >= start & pos <= end, 0,
                                    pmin(abs(pos - start), abs(pos - end)))) %>%
        dplyr::arrange(dist)
      return(c(bin = hits$bin[1], gene = hits$gene[1], how = "cis-window"))
    }
    if (!is.na(ng) && ng %in% panel_sym$gene) {
      return(c(bin = panel_sym$bin[match(ng, panel_sym$gene)], gene = ng, how = "nearest-gene"))
    }
    c(bin = "other_pleiotropic", gene = NA, how = "none")
  }
  ng_vec <- if (gene_col %in% names(df)) df[[gene_col]] else rep(NA_character_, nrow(df))
  res <- purrr::pmap(list(df[[chr_col]], df[[pos_col]], ng_vec), assign_one)
  df %>% mutate(mechanism   = vapply(res, `[[`, character(1), "bin"),
                bin_gene    = vapply(res, `[[`, character(1), "gene"),
                bin_method  = vapply(res, `[[`, character(1), "how"))
}

# ---- Instrument strength -----------------------------------------------------
# Per-SNP R2 and F, plus an overall summary tibble. Expects harmonised TwoSampleMR cols.
add_instrument_strength <- function(df, n_exposure) {
  df %>% mutate(
    samplesize.exposure = n_exposure,
    R2.exposure = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
    F.exposure  = (R2.exposure * (samplesize.exposure - 2)) / (1 - R2.exposure)
  )
}

instrument_summary <- function(df) {
  df %>% summarise(
    num_snps = n(),
    samplesize.exposure = dplyr::first(samplesize.exposure),
    cumulative_R2 = sum(R2.exposure, na.rm = TRUE),
    mean_F = mean(F.exposure, na.rm = TRUE),
    median_F = median(F.exposure, na.rm = TRUE)
  ) %>% mutate(
    overall_F = (cumulative_R2 * (samplesize.exposure - num_snps - 1)) /
                ((1 - cumulative_R2) * num_snps)
  )
}

# ---- Relaxed cis drug-target MR (LD-aware) -----------------------------------
# Select relaxed cis instruments for an exposure within a physical window, allowing
# correlated SNPs (clump_r2), then run LD-aware MR (IVW + Egger) against an outcome using
# the variant correlation matrix. Convention follows cholesterol drug-target cis-MR.

# Retry an OpenGWAS call with exponential backoff (handles transient rate-limiting).
# `fn` must be a zero-arg function (thunk) so it can be re-evaluated on each attempt.
# Returns the result, or NULL if it never yields a non-empty data frame.
og_retry <- function(fn, tries = 4, base_sleep = 3) {
  res <- NULL
  for (i in seq_len(tries)) {
    res <- tryCatch(fn(), error = function(e) NULL)
    ok <- !is.null(res) && (!is.data.frame(res) || nrow(res) > 0)
    if (ok) return(res)
    Sys.sleep(base_sleep * i)
  }
  res
}

# Pull + relax-clump cis instruments for `exposure_id` in `region` ("chr:start-end").
# Returns a TwoSampleMR-formatted exposure tibble, or NULL if none.
cis_instruments_relaxed <- function(exposure_id, region, cfg) {
  suppressPackageStartupMessages(library(ieugwasr))
  a <- og_retry(function() associations(region, exposure_id))
  if (is.null(a) || nrow(a) == 0) return(NULL)
  a <- a |> dplyr::filter(p < cfg$cis_drugtarget$p_thresh)
  if (nrow(a) == 0) return(NULL)
  if (nrow(a) > 1) {
    cl <- og_retry(function()
      ld_clump(dplyr::tibble(rsid = a$rsid, pval = a$p, id = a$id),
               clump_r2 = cfg$cis_drugtarget$clump_r2,
               clump_kb = cfg$cis_drugtarget$clump_kb,
               pop = cfg$cis_drugtarget$ld_pop))
    if (!is.null(cl)) a <- a |> dplyr::filter(rsid %in% cl$rsid)
  }
  a |> dplyr::transmute(
    SNP = rsid, beta.exposure = beta, se.exposure = se,
    effect_allele.exposure = ea, other_allele.exposure = nea,
    eaf.exposure = eaf, pval.exposure = p,
    exposure = exposure_id, id.exposure = exposure_id,
    mr_keep.exposure = TRUE, pval_origin.exposure = "reported",
    chr.exposure = chr, pos.exposure = position)
}

# LD-aware cis-MR of a formatted exposure set against `outcome_id`.
# Harmonises, fetches the LD matrix, aligns alleles, runs MendelianRandomization IVW+Egger
# with the correlation matrix. Returns a tidy tibble (one row per method) + Q heterogeneity.
ld_aware_cis_mr <- function(exp_df, outcome_id, label, cfg) {
  suppressPackageStartupMessages({
    library(TwoSampleMR); library(ieugwasr); library(MendelianRandomization)
  })
  empty <- tibble(locus = label, method = NA_character_, nsnp = 0L,
                  b = NA_real_, se = NA_real_, pval = NA_real_)
  if (is.null(exp_df) || nrow(exp_df) == 0) return(empty)
  out <- og_retry(function() extract_outcome_data(exp_df$SNP, outcome_id, proxies = TRUE, rsq = 0.8))
  if (is.null(out) || nrow(out) == 0) return(empty)
  h <- suppressMessages(harmonise_data(exp_df, out, action = 2)) |> dplyr::filter(mr_keep)
  if (nrow(h) == 0) return(empty)
  if (nrow(h) == 1) {                       # single SNP -> Wald ratio, no LD needed
    b <- h$beta.outcome / h$beta.exposure
    se <- abs(h$se.outcome / h$beta.exposure)
    return(tibble(locus = label, method = "Wald ratio", nsnp = 1L,
                  b = b, se = se, pval = 2 * pnorm(-abs(b / se))))
  }
  ld <- og_retry(function() ld_matrix(h$SNP, with_alleles = TRUE, pop = cfg$cis_drugtarget$ld_pop))
  if (is.null(ld)) return(empty)
  hl <- tryCatch(TwoSampleMR::harmonise_ld_dat(h, ld), error = function(e) NULL)
  if (is.null(hl) || nrow(hl$x) < 2) return(empty)
  dat <- hl$x; ldm <- hl$ld
  mri <- MendelianRandomization::mr_input(
    bx = dat$beta.exposure, bxse = dat$se.exposure,
    by = dat$beta.outcome,  byse = dat$se.outcome, correlation = ldm)
  ivw <- tryCatch(MendelianRandomization::mr_ivw(mri, correl = TRUE), error = function(e) NULL)
  egg <- tryCatch(MendelianRandomization::mr_egger(mri, correl = TRUE), error = function(e) NULL)
  res <- list()
  if (!is.null(ivw)) res[[length(res)+1]] <- tibble(
    locus = label, method = "LD-aware IVW", nsnp = ivw$SNPs,
    b = ivw$Estimate, se = ivw$StdError, pval = ivw$Pvalue)
  if (!is.null(egg)) res[[length(res)+1]] <- tibble(
    locus = label, method = "LD-aware MR-Egger", nsnp = egg$SNPs,
    b = egg$Estimate, se = egg$StdError.Est, pval = egg$Pvalue.Est)
  if (length(res) == 0) return(empty)
  dplyr::bind_rows(res)
}

# ---- Equivalence testing (TOST) for selection negative controls --------------
# Two one-sided tests of equivalence against a symmetric bound +/- delta on the MR
# estimate (b, se). Returns the TOST p-value (max of the two one-sided p's) plus the
# 95% CI and the pre-registered verdict. Per APPROACH.md sec 13 a cell is "equiv-null"
# only if the 95% CI is fully inside +/- delta (stricter than TOST p<0.05, which maps
# to the 90% CI); a significant association (p_assoc<0.05) not inside the bound is
# "real-assoc"; anything else is "inconclusive".
tost_equivalence <- function(b, se, delta) {
  if (is.na(b) || is.na(se) || se <= 0)
    return(tibble::tibble(tost_p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                          p_assoc = NA_real_, equiv_ci95 = NA, verdict = NA_character_))
  t_l <- (b + delta) / se; p_l <- stats::pnorm(t_l, lower.tail = FALSE)  # H0: b <= -delta
  t_u <- (b - delta) / se; p_u <- stats::pnorm(t_u, lower.tail = TRUE)   # H0: b >= +delta
  tost_p <- max(p_l, p_u)
  ci_lo <- b - 1.96 * se; ci_hi <- b + 1.96 * se
  p_assoc <- 2 * stats::pnorm(-abs(b / se))
  equiv_ci95 <- (ci_lo >= -delta) && (ci_hi <= delta)
  verdict <- if (isTRUE(equiv_ci95)) "equiv-null"
             else if (!is.na(p_assoc) && p_assoc < 0.05) "real-assoc"
             else "inconclusive"
  tibble::tibble(tost_p = tost_p, ci_lo = ci_lo, ci_hi = ci_hi,
                 p_assoc = p_assoc, equiv_ci95 = equiv_ci95, verdict = verdict)
}

# Pull + cache the relaxed cis instrument set for the smoking exposure at one nAChR locus.
# Caches the formatted exposure tibble under data/cache/instr_<locus>.rds so Task A/B reuse
# the SAME fixed instruments without re-hitting OpenGWAS. Returns the exposure tibble (or NULL).
get_locus_instruments <- function(locus_name, cfg, refresh = FALSE) {
  cache_dir <- here::here("data", "cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  f <- file.path(cache_dir, paste0("instr_", locus_name, ".rds"))
  if (!refresh && file.exists(f)) return(readRDS(f))
  l <- cfg$nachr_loci[[locus_name]]
  win_kb <- cfg$cis_window_kb
  region <- sprintf("%d:%d-%d", l$chr, l$start - win_kb * 1000, l$end + win_kb * 1000)
  exp <- cis_instruments_relaxed(cfg$opengwas$smk_cpd, region, cfg)
  if (!is.null(exp)) saveRDS(exp, f)
  exp
}

# ---- Compile-safe results: never clobber a good cache with a failed live pull ---------------
# If the live result has at least one non-NA value in `key_col`, write it to `csv` and return
# it. Otherwise (OpenGWAS unreachable / token expired at compile time) fall back to the cached
# csv if present, so a failed pull cannot overwrite verified results and downstream plots/tables
# still render. Pass the live result (or NULL) as `live`.
# One-shot, session-memoized OpenGWAS connectivity probe. Lets a chunk skip the (slow,
# back-off-retrying) live pull entirely when the token is expired/quota-exhausted at compile
# time, so the doc falls straight back to cached results instead of hanging for minutes.
.opengwas_ok_cache <- new.env(parent = emptyenv())
opengwas_ok <- function(refresh = FALSE) {
  if (!refresh && !is.null(.opengwas_ok_cache$ok)) return(.opengwas_ok_cache$ok)
  ok <- tryCatch(nrow(ieugwasr::gwasinfo("ieu-b-142")) > 0, error = function(e) FALSE)
  if (!ok) message("OpenGWAS not reachable (token/quota) — layers will use cached results.")
  .opengwas_ok_cache$ok <- ok
  ok
}

live_or_cache <- function(live, csv, key_col) {
  good <- !is.null(live) && (key_col %in% names(live)) && any(!is.na(live[[key_col]]))
  if (good) { readr::write_csv(live, csv); return(live) }
  if (file.exists(csv)) {
    message("Live pull empty/failed — using cached ", basename(csv), " (not overwritten)")
    return(readr::read_csv(csv, show_col_types = FALSE))
  }
  live
}

# ---- Decisions log -----------------------------------------------------------
# Append a single decision line to results/decisions_log.md.
log_decision <- function(layer, rule, value, decision,
                         path = here::here("results", "decisions_log.md")) {
  if (!file.exists(path)) {
    writeLines(c("# Decisions log",
                 "",
                 "| Layer | Rule | Value | Decision |",
                 "|---|---|---|---|"), path)
  }
  cat(sprintf("| %s | %s | %s | %s |\n", layer, rule, value, decision),
      file = path, append = TRUE)
}
