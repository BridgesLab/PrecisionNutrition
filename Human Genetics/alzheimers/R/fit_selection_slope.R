# fit_selection_slope.R — generalized SlopeHunter selection-slope engine.
#
# Estimates the survival/selection-induced slope b_SH that the AD outcome carries because it is
# measured only in survivors, and provides the exact linear correction that removes it from any
# exposure->AD estimate. This is the multi-exposure generalization of the smoking pipeline's
# smoking/R/prep_indexevent_genomewide.R + run_indexevent_genomewide.R, with the smoking-specific
# 15q25 handling dropped (there is no single "locus under test" across a battery of exposures).
#
# The selection axis (SlopeHunter's "incidence") is a lifespan/survival GWAS; the "prognosis" is
# AD. b_SH relates SNP effects on AD to SNP effects on lifespan among survivors, so it is a
# property of the (AD GWAS x lifespan GWAS) pair and is exposure-independent: estimate it once
# per AD GWAS, then apply to every exposure.
#
# Because IVW is linear, the corrected exposure->AD slope is exact:
#   corrected(E->AD) = (E->AD) - b_SH * (E->lifespan)
#
# Heavy inputs: this needs FULL summary statistics for AD and lifespan (the OpenGWAS API only
# serves clumped tophits). Reuse smoking/scripts/fetch_sumstats.sh to download them, and plink2
# --clump for the LD-independent SNP set. The .qmd consumes the fits CSV this writes; it does not
# refit inline.

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})

# ---- allele helpers ----------------------------------------------------------
is_palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
is_snp <- function(a1, a2) nchar(a1) == 1L & nchar(a2) == 1L &
  a1 %in% c("A", "C", "G", "T") & a2 %in% c("A", "C", "G", "T")

# ---- Step 1: build the allele-aligned AD x lifespan merge --------------------
# ad_file / life_file: gzipped GWAS-Catalog-harmonised TSVs (variant_id, chromosome,
#   base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta,
#   standard_error, p_value). AD effects are aligned to the LIFESPAN effect allele (SlopeHunter
#   treats lifespan as x). Palindromic / multiallelic / allele-mismatch SNPs are dropped.
prep_selection_merge <- function(ad_file, life_file, life_p_keep = 0.01) {
  cols <- c("variant_id", "chromosome", "base_pair_location", "effect_allele",
            "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")

  message("Reading lifespan (incidence/selection axis) ...")
  life <- fread(cmd = paste("gunzip -c", shQuote(life_file)), select = cols, showProgress = FALSE)
  setnames(life, c("SNP","chr_life","pos_life","ea_life","oa_life","eaf_life",
                   "beta_life","se_life","p_life"))
  life <- life[!is.na(p_life) & p_life < life_p_keep &
                 !is.na(beta_life) & !is.na(se_life) & se_life > 0]
  life <- life[is_snp(ea_life, oa_life)]
  message("  lifespan SNPs with p < ", life_p_keep, ": ", nrow(life))

  message("Reading AD (prognosis) ...")
  ad <- fread(cmd = paste("gunzip -c", shQuote(ad_file)), select = cols, showProgress = FALSE)
  setnames(ad, c("SNP","chr_ad","pos_ad","ea_ad","oa_ad","eaf_ad","beta_ad","se_ad","p_ad"))
  ad <- ad[!is.na(beta_ad) & !is.na(se_ad) & se_ad > 0][is_snp(ea_ad, oa_ad)]
  message("  AD SNPs usable: ", nrow(ad))

  m <- merge(life, ad, by = "SNP")
  m <- m[!is_palindromic(ea_life, oa_life)]
  same <- m$ea_ad == m$ea_life & m$oa_ad == m$oa_life
  flip <- m$ea_ad == m$oa_life & m$oa_ad == m$ea_life
  m <- m[same | flip]
  m[, beta_ad_aligned := fifelse(ea_ad == ea_life, beta_ad, -beta_ad)]
  m <- m[is.na(eaf_life) | (eaf_life > 0.01 & eaf_life < 0.99)]

  # APOE flag in both builds (AD often GRCh38, lifespan often GRCh37) — genuine AD-longevity
  # pleiotropy, excluded from the PRIMARY fit because it is not a selection collider.
  m[, apoe := (chr_ad == 19 & pos_ad > 43.9e6 & pos_ad < 45.9e6) |
              (chr_life == 19 & pos_life > 44.4e6 & pos_life < 46.5e6)]
  message("  merged/aligned: ", nrow(m), " | APOE-region: ", sum(m$apoe))
  m[]
}

# ---- Step 2: fit SlopeHunter on one SNP subset -------------------------------
fit_slopehunter_one <- function(d, label, xp_thresh = 1e-3, seed = 20260604) {
  suppressPackageStartupMessages(library(SlopeHunter))
  dat <- data.frame(
    SNP            = d$SNP,
    BETA.incidence = d$beta_life,        SE.incidence = d$se_life, Pval.incidence = d$p_life,
    BETA.prognosis = d$beta_ad_aligned,  SE.prognosis = d$se_ad,   Pval.prognosis = d$p_ad,
    stringsAsFactors = FALSE)
  n_fit <- sum(dat$Pval.incidence < xp_thresh, na.rm = TRUE)
  fit <- tryCatch(SlopeHunter::hunt(dat, xp_thresh = xp_thresh, Plot = FALSE, seed = seed),
                  error = function(e) { message("hunt() failed for ", label, ": ",
                                                conditionMessage(e)); NULL })
  if (is.null(fit))
    return(tibble(set = label, n_snp = nrow(dat), n_fit = n_fit,
                  b_SH = NA_real_, se_SH = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  b  <- as.numeric(fit$b)
  se <- suppressWarnings(as.numeric(fit$bse))
  if (length(se) == 0 || is.na(se)) se <- stats::sd(fit$Bts.est, na.rm = TRUE)
  tibble(set = label, n_snp = nrow(dat), n_fit = n_fit,
         b_SH = b, se_SH = se, ci_lo = b - 1.96*se, ci_hi = b + 1.96*se)
}

# Fit the full set and the APOE-excluded PRIMARY set.
fit_selection_slope <- function(m, xp_thresh = 1e-3, seed = 20260604) {
  bind_rows(
    fit_slopehunter_one(m,                "full (incl APOE)",  xp_thresh, seed),
    fit_slopehunter_one(m[apoe == FALSE], "excl_APOE (PRIMARY)", xp_thresh, seed))
}

# ---- exact linear correction (delta-method SE) -------------------------------
# corrected(E->AD) = b_ad - b_SH * b_life. Reused by the .qmd.
correct_estimate <- function(b_ad, se_ad, b_life, se_life, b_SH, se_SH) {
  b_corr  <- b_ad - b_SH * b_life
  se_corr <- sqrt(se_ad^2 + (b_SH^2) * se_life^2 + (b_life^2) * se_SH^2)
  tibble(b_uncorr = b_ad, b_corr = b_corr, se_corr = se_corr,
         ci_lo = b_corr - 1.96*se_corr, ci_hi = b_corr + 1.96*se_corr,
         p_corr = 2 * stats::pnorm(-abs(b_corr / se_corr)))
}

# ---- orientation diagnostic ------------------------------------------------
# Aligned AD betas for the SAME disease should correlate POSITIVELY across two AD GWAS (both are
# aligned to the lifespan effect allele in prep_selection_merge). A NEGATIVE correlation means one
# GWAS's effect alleles are systematically inverted vs the other, which would flip its b_SH sign
# (and any correction built on it). Use this when a new AD GWAS's b_SH comes out opposite-signed
# to a trusted one (e.g. Kunkle +1.10 vs Bellenguez -0.944). Pass two prep_selection_merge()
# outputs; the first is the trusted reference.
ad_orientation_check <- function(m_ref, m_test, ref_label = "ref", test_label = "test") {
  j <- merge(m_ref[,  .(SNP, beta_ref  = beta_ad_aligned)],
             m_test[, .(SNP, beta_test = beta_ad_aligned)], by = "SNP")
  if (nrow(j) < 20)
    return(tibble::tibble(ref = ref_label, test = test_label, n_shared = nrow(j),
                          correlation = NA_real_, sign_agreement = NA_real_,
                          verdict = "too few shared SNPs"))
  r     <- stats::cor(j$beta_ref, j$beta_test)
  agree <- mean(sign(j$beta_ref) == sign(j$beta_test))
  tibble::tibble(
    ref = ref_label, test = test_label, n_shared = nrow(j),
    correlation = r, sign_agreement = agree,
    verdict = if (r < 0)                     "LIKELY FLIPPED — test AD alleles inverted vs ref"
              else if (r > 0.5 & agree > 0.7) "consistent orientation"
              else                            "ambiguous — inspect")
}

# ---- orientation diagnostic 2: single-file anchor check ----------------------
# ad_orientation_check() needs TWO prepped merges, so it cannot run when only one AD sumstats
# file is on disk. This version needs just the one file: it reads published AD risk alleles at
# well-established loci and asks whether the file's betas agree. APOE rs429358-C alone is
# decisive (p ~ 1e-881 in Kunkle) — no correctly-oriented AD GWAS can get its sign wrong.
# Returns one row per anchor plus an overall verdict attribute.
AD_ANCHORS <- tibble::tribble(
  ~SNP,         ~risk_allele, ~locus,
  "rs429358",   "C",          "APOE e4 (decisive)",
  "rs7412",     "C",          "APOE (non-e2)",
  "rs6656401",  "A",          "CR1",
  "rs6733839",  "T",          "BIN1",
  "rs11136000", "C",          "CLU",
  "rs3851179",  "C",          "PICALM",
  "rs10948363", "G",          "CD2AP"
)

check_ad_orientation_anchors <- function(ad_file, anchors = AD_ANCHORS) {
  cols <- c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "p_value")
  ad <- data.table::fread(cmd = paste("gunzip -c", shQuote(ad_file)), select = cols,
                          showProgress = FALSE)
  data.table::setnames(ad, c("SNP", "ea", "oa", "beta", "se", "p"))
  ad <- ad[SNP %in% anchors$SNP]
  res <- dplyr::inner_join(anchors, tibble::as_tibble(ad), by = "SNP") |>
    dplyr::mutate(
      # beta re-expressed for the PUBLISHED risk allele; must be > 0 in a correct file.
      beta_risk = dplyr::if_else(ea == risk_allele, beta, -beta),
      ok        = beta_risk > 0)
  n_ok <- sum(res$ok, na.rm = TRUE)
  apoe <- res |> dplyr::filter(SNP == "rs429358")
  verdict <-
    if (nrow(apoe) == 1 && !apoe$ok[1]) "FLIPPED — rs429358 risk allele has a negative beta"
    else if (n_ok == nrow(res))         "correct orientation (all anchors agree)"
    else if (n_ok >= nrow(res) - 1)     "correct orientation (APOE + majority agree)"
    else                                "ambiguous — inspect"
  structure(res, verdict = verdict)
}

# ---- b_SH identifiability diagnostic ----------------------------------------
# SlopeHunter ALWAYS returns a slope, with a bootstrap SE that looks precise, even when the two
# axes are unrelated: its mixture model partitions whatever cloud it is given, and with no signal
# the "hunted" cluster is carved out of noise. So a b_SH is only interpretable if the underlying
# data actually carry a selection relationship. This reports the raw evidence for one, on the
# SAME SNPs SlopeHunter fits (incidence p < xp_thresh). Near-zero r with ~0.5 sign agreement
# means b_SH is not estimable from these data regardless of how tight its CI looks.
bsh_identifiability <- function(m, xp_thresh = 1e-3) {
  d <- m[p_life < xp_thresh]
  r  <- stats::cor(d$beta_life, d$beta_ad_aligned)
  ols <- stats::coef(stats::lm(beta_ad_aligned ~ 0 + beta_life, data = d))[[1]]
  tibble::tibble(
    n_fit          = nrow(d),
    pearson_r      = r,
    ols_slope      = ols,
    sign_agreement = mean(sign(d$beta_life) == sign(d$beta_ad_aligned)),
    sd_ratio       = stats::sd(d$beta_ad_aligned) / stats::sd(d$beta_life),
    verdict        = if (abs(r) < 0.05)
                       "NOT IDENTIFIED — no raw selection signal; b_SH reflects noise geometry"
                     else if (abs(r) < 0.15) "weak — treat b_SH as provisional"
                     else                    "signal present — b_SH interpretable")
}

# ---- null distribution of b_SH ----------------------------------------------
# The definitive check on a b_SH. bsh_identifiability() shows whether raw signal exists; this
# shows what SlopeHunter returns when it provably does NOT. The AD block (beta, se, p) is
# permuted TOGETHER across SNPs, so each SNP keeps internally consistent AD statistics while the
# lifespan<->AD pairing is destroyed. Whatever comes back is pure artefact of the data geometry.
#
# The key result: THE NULL IS NOT CENTRED ON ZERO. On the Kunkle merge, permuted data return
# b_SH ~ -0.73 (range -0.73 to -0.53, sd 0.074) — a large negative slope produced by nothing but
# the noise geometry of the two axes. The bootstrap SE (0.112) is not misleadingly narrow; it is
# the *reference point* that is wrong. So "b_SH differs significantly from 0" is NOT evidence of
# selection. Compare an observed b_SH against this permutation null, never against zero.
#
# Slow (~10-20 s per replicate); n_perm = 8 is enough to see the spread.
bsh_null_permutation <- function(m, n_perm = 8, xp_thresh = 1e-3, seed = 42) {
  suppressPackageStartupMessages(library(SlopeHunter))
  set.seed(seed)
  b <- vapply(seq_len(n_perm), function(i) {
    k <- sample(nrow(m))
    dat <- data.frame(
      SNP            = m$SNP,
      BETA.incidence = m$beta_life,           SE.incidence = m$se_life,
      Pval.incidence = m$p_life,
      BETA.prognosis = m$beta_ad_aligned[k],  SE.prognosis = m$se_ad[k],
      Pval.prognosis = m$p_ad[k], stringsAsFactors = FALSE)
    f <- tryCatch(SlopeHunter::hunt(dat, xp_thresh = xp_thresh, Plot = FALSE, seed = 7),
                  error = function(e) NULL)
    if (is.null(f)) NA_real_ else as.numeric(f$b)
  }, numeric(1))
  tibble::tibble(n_perm = n_perm, null_min = min(b, na.rm = TRUE),
                 null_max = max(b, na.rm = TRUE), null_median = stats::median(b, na.rm = TRUE),
                 null_sd = stats::sd(b, na.rm = TRUE))
}

# ---- driver: produce a fits CSV for one AD GWAS ------------------------------
# clumped_snps: optional character vector (or path to a one-column file) of LD-independent rsIDs
#   from plink2 --clump. If NULL, fits on all merged SNPs (not recommended — LD inflates n).
build_selection_slope <- function(ad_file, life_file, out_csv,
                                   clumped_snps = NULL, life_p_keep = 0.01,
                                   xp_thresh = 1e-3, seed = 20260604) {
  m <- prep_selection_merge(ad_file, life_file, life_p_keep = life_p_keep)
  if (is.character(clumped_snps) && length(clumped_snps) == 1 && file.exists(clumped_snps)) {
    ct <- data.table::fread(clumped_snps, header = TRUE)
    idcol <- intersect(c("ID", "SNP", "RSID", "rsid", "MarkerName"), names(ct))
    clumped_snps <- if (length(idcol)) ct[[idcol[1]]] else ct[[1]]  # plink2 .clumps -> "ID"
  }
  if (!is.null(clumped_snps)) {
    m <- m[SNP %in% clumped_snps]
    message("SNPs after LD clumping: ", nrow(m))
  }
  fits <- fit_selection_slope(m, xp_thresh = xp_thresh, seed = seed)
  if (!is.null(out_csv)) readr::write_csv(fits, out_csv)
  fits
}
