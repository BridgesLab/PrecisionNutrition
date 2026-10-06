# read.R — read full GWAS summary statistics into one standard format, with a QC attrition log.
#
# Standard format (data.table, one row per variant, effect allele = ea):
#   SNP chr pos ea oa eaf beta se p n
# Every downstream step works on this, so method differences are never input-format differences.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

`%||%` <- function(a, b) if (is.null(a)) b else a

# OpenGWAS GWAS-VCF. ES/SE/LP/AF/SS are FORMAT fields and the FORMAT string varies by row (AF and
# SS are present on some rows only), so rows are split per distinct FORMAT string. ALT is the
# effect allele; LP is -log10(p).
#
# `keep_snps` (the LD panel's rsIDs) is applied while streaming, in awk, before R sees the file:
# the 2hGlu VCF has 27M rows, which would need ~19 GB in R. Variants outside the panel cannot be
# clumped or given an LD matrix, so no method could use them anyway. The total row count is kept
# (attr "n_total") so the attrition table still starts from the full file.
read_gwas_vcf <- function(path, keep_snps = NULL) {
  stream <- paste("gunzip -c", shQuote(path), "| grep -v '^##'")
  n_total <- NA_real_
  if (!is.null(keep_snps)) {
    keyf <- tempfile(fileext = ".txt"); cnt <- tempfile(); on.exit(unlink(c(keyf, cnt)))
    writeLines(keep_snps, keyf)
    stream <- paste0(stream, " | awk -F'\\t' -v cnt=", shQuote(cnt),
                     " 'NR==FNR{k[$1];next} /^#/{print;next} {n++} ($3 in k){print}",
                     " END{print n > cnt}' ", shQuote(keyf), " -")
  }
  d <- fread(cmd = stream, sep = "\t", header = TRUE, select = c(1, 2, 3, 4, 5, 9, 10),
             colClasses = "character", showProgress = FALSE)
  if (!is.null(keep_snps)) n_total <- as.numeric(readLines(cnt))
  setnames(d, c("chr", "pos", "SNP", "oa", "ea", "fmt", "val"))
  out <- d[, {
    keys <- strsplit(fmt[1], ":", fixed = TRUE)[[1]]
    v <- tstrsplit(val, ":", fixed = TRUE)
    names(v) <- keys
    g <- function(k) if (k %in% keys) suppressWarnings(as.numeric(v[[k]])) else NA_real_
    list(SNP = SNP, chr = chr, pos = pos, ea = ea, oa = oa,
         beta = g("ES"), se = g("SE"), lp = g("LP"), eaf = g("AF"), n = g("SS"))
  }, by = fmt][, fmt := NULL]
  out[, `:=`(p = 10^(-lp), pos = as.integer(pos))][, lp := NULL]
  setattr(out, "n_total", n_total)
  out[]
}

# Local GWAS-Catalog-harmonised TSV (the Kunkle and lifespan files Step 2 already uses).
read_gwas_tsv <- function(path) {
  cols <- c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "effect_allele_frequency", "beta", "standard_error", "p_value")
  d <- fread(cmd = paste("gunzip -c", shQuote(path)), select = cols, showProgress = FALSE)
  setnames(d, c("SNP", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p"))
  # p-values below the double range (Kunkle APOE: 1.2e-881) make fread read the column as text;
  # as.numeric() turns those into 0, which qc_sumstats() floors rather than drops.
  d[, `:=`(chr = as.character(chr), n = NA_real_, eaf = as.numeric(eaf),
           p = suppressWarnings(as.numeric(p)))]
  d[]
}

# QC to the analysis set. Each step is counted so the attrition table can show where variants went.
#   - rsID present; biallelic SNV (single A/C/G/T alleles)
#   - finite beta, se > 0, p in (0, 1]
#   - MAF >= maf_min where a frequency is available (missing frequency is kept, as OpenGWAS did)
#   - one row per rsID (lowest p kept)
#   - present in the LD reference panel (needed by every clumping / LD step)
# `n_fixed` fills a missing per-SNP sample size; `p_from_beta` recomputes p from beta/se, for files
# whose p-value column comes from a different model than the reported beta (2hGlu: INT p, raw beta).
qc_sumstats <- function(d, dataset, ref_snps, maf_min = 0.01, n_fixed = NA_real_, p_from_beta = FALSE) {
  steps <- list()
  note <- function(step) steps[[length(steps) + 1]] <<- tibble(dataset = dataset, step = step,
                                                                n_variants = nrow(d))
  if (!is.null(attr(d, "n_total")) && !is.na(attr(d, "n_total")))
    steps[[1]] <- tibble(dataset = dataset, step = "variants in file", n_variants = attr(d, "n_total"))
  note(if (length(steps)) "rsID in LD reference panel (filtered while reading)" else "read")
  d <- d[grepl("^rs[0-9]+$", SNP)];                                     note("rsID present")
  d[, `:=`(ea = toupper(ea), oa = toupper(oa))]
  d <- d[ea %chin% c("A", "C", "G", "T") & oa %chin% c("A", "C", "G", "T")]; note("biallelic SNV")
  if (p_from_beta) d[, p := 2 * pnorm(-abs(beta / se))]
  # Underflowed p (0) is real and extreme, not missing: floor it so the variant is kept.
  n_floor <- d[is.finite(p) & p == 0, .N]
  d[is.finite(p) & p == 0, p := .Machine$double.xmin]
  d <- d[is.finite(beta) & is.finite(se) & se > 0 & is.finite(p) & p > 0 & p <= 1]
  note(paste0(if (p_from_beta) "valid beta/se (p recomputed from beta/se)" else "valid beta/se/p",
              if (n_floor) sprintf(" [%d p-values below double range floored, not dropped]", n_floor) else ""))
  d <- d[is.na(eaf) | (eaf >= maf_min & eaf <= 1 - maf_min)];           note(paste0("MAF >= ", maf_min))
  setorder(d, p)
  d <- unique(d, by = "SNP");                                           note("unique rsID")
  d <- d[SNP %chin% ref_snps];                                          note("in LD reference panel")
  if (!is.na(n_fixed)) d[is.na(n), n := n_fixed]
  setattr(d, "qc_log", bind_rows(steps))
  d[]
}

# Dispatcher used by the pipeline: source type from config -> standard format -> QC.
load_gwas <- function(path, ds, ref_snps, qc_cfg) {
  raw <- switch(ds$source,
                opengwas_vcf = read_gwas_vcf(path, keep_snps = ref_snps),
                local_tsv    = read_gwas_tsv(path),
                stop("unknown source '", ds$source, "' for ", ds$key))
  qc_sumstats(raw, ds$key, ref_snps, maf_min = qc_cfg$maf_min, n_fixed = ds$n %||% NA_real_,
          p_from_beta = isTRUE(ds$p_from_beta))
}

# SNP set of the plink panel (bim column 2) with its alleles; col5 is the allele plink2 counts.
read_bim <- function(bfile) {
  b <- fread(paste0(bfile, ".bim"), header = FALSE, select = c(1, 2, 4, 5, 6),
             col.names = c("chr", "SNP", "pos", "a1", "a2"), showProgress = FALSE)
  b[!duplicated(SNP)]
}
