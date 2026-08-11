# run_indexevent_genomewide.R — Step 2 of the properly-powered index-event correction.
#
# Fits SlopeHunter on the genome-wide, LD-clumped, allele-aligned merge produced by
# prep_indexevent_genomewide.R, then re-applies the linear correction to the nAChR estimates.
#
# Three fits are reported so the contamination issues that broke the original n=10 run are
# visible rather than buried:
#   full            - all clumped SNPs (comparable to the old, contaminated run)
#   excl_APOE       - PRIMARY: APOE removed (genuine AD-longevity pleiotropy, not a collider)
#   excl_APOE_15q25 - SENSITIVITY: also removes the locus under test (circularity)
#
# Because IVW is linear, the corrected smoking->AD slope is exact:
#   corrected(smk->AD) = (smk->AD) - b_SH * (smk->lifespan)

suppressPackageStartupMessages({
  library(data.table); library(SlopeHunter); library(here); library(dplyr)
})

# Fit SlopeHunter on one SNP subset. Returns slope, SE and n.
fit_slopehunter <- function(d, label, xp_thresh = 1e-3, seed = 20260604) {
  dat <- data.frame(
    SNP            = d$SNP,
    BETA.incidence = d$beta_life,  SE.incidence = d$se_life,  Pval.incidence = d$p_life,
    BETA.prognosis = d$beta_ad_aligned, SE.prognosis = d$se_ad, Pval.prognosis = d$p_ad,
    stringsAsFactors = FALSE)
  n_fit <- sum(dat$Pval.incidence < xp_thresh, na.rm = TRUE)
  fit <- tryCatch(
    SlopeHunter::hunt(dat, xp_thresh = xp_thresh, Plot = FALSE, seed = seed),
    error = function(e) { message("hunt() failed for ", label, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(tibble::tibble(set = label, n_snp = nrow(dat), n_fit = n_fit,
                                          b_SH = NA_real_, se_SH = NA_real_,
                                          ci_lo = NA_real_, ci_hi = NA_real_))
  b <- as.numeric(fit$b)
  se <- suppressWarnings(as.numeric(fit$bse))
  if (length(se) == 0 || is.na(se)) se <- stats::sd(fit$Bts.est, na.rm = TRUE)
  tibble::tibble(set = label, n_snp = nrow(dat), n_fit = n_fit,
                 b_SH = b, se_SH = se, ci_lo = b - 1.96*se, ci_hi = b + 1.96*se)
}

# Apply the exact linear correction to a smoking->AD IVW estimate.
correct_estimate <- function(b_ad, se_ad, b_life, se_life, b_SH, se_SH) {
  b_corr  <- b_ad - b_SH * b_life
  se_corr <- sqrt(se_ad^2 + (b_SH^2) * se_life^2 + (b_life^2) * se_SH^2)
  tibble::tibble(b_uncorr = b_ad, b_corr = b_corr, se_corr = se_corr,
                 ci_lo = b_corr - 1.96*se_corr, ci_hi = b_corr + 1.96*se_corr,
                 p_corr = 2*stats::pnorm(-abs(b_corr/se_corr)))
}

run_indexevent_genomewide <- function(clumped_snps = NULL) {
  m <- readRDS(here::here("data","cache","sumstats","merged_AD_lifespan.rds"))
  if (!is.null(clumped_snps)) m <- m[SNP %in% clumped_snps]
  message("SNPs entering SlopeHunter (clumped): ", nrow(m))

  fits <- dplyr::bind_rows(
    fit_slopehunter(m,                                    "full (incl APOE + 15q25)"),
    fit_slopehunter(m[apoe == FALSE],                     "excl_APOE (PRIMARY)"),
    fit_slopehunter(m[apoe == FALSE & locus15q25 == FALSE], "excl_APOE_15q25 (sensitivity)"))

  # nAChR smoking->AD (L6) and smoking->lifespan (Task A) inputs for the correction.
  pooled_ad <- readr::read_csv(here::here("results","cis_smoking_ldaware_pooled.csv"),
                               show_col_types = FALSE) |> dplyr::filter(method == "LD-aware IVW")
  per15 <- readr::read_csv(here::here("results","cis_smoking_ldaware_per_locus.csv"),
                           show_col_types = FALSE) |>
    dplyr::filter(locus == "CHRNA5_A3_B4", method == "LD-aware IVW")
  sb <- readr::read_csv(here::here("results","selection_battery.csv"), show_col_types = FALSE) |>
    dplyr::filter(outcome == "parental_lifespan")

  targets <- list(
    list(set = "pooled_nAChR",     b = pooled_ad$b, se = pooled_ad$se,
         bl = sb$b[sb$set == "pooled_nAChR"],     sl = sb$se[sb$set == "pooled_nAChR"]),
    list(set = "15q25_CHRNA5A3B4", b = per15$b,     se = per15$se,
         bl = sb$b[sb$set == "15q25_CHRNA5A3B4"], sl = sb$se[sb$set == "15q25_CHRNA5A3B4"]))

  corr <- purrr::map_dfr(seq_len(nrow(fits)), function(i) {
    f <- fits[i, ]
    purrr::map_dfr(targets, function(tg) {
      if (is.na(f$b_SH)) return(NULL)
      correct_estimate(tg$b, tg$se, tg$bl, tg$sl, f$b_SH, f$se_SH) |>
        dplyr::mutate(slopehunter_set = f$set, nachr_set = tg$set,
                      b_SH = f$b_SH, se_SH = f$se_SH, n_fit = f$n_fit, .before = 1)
    })
  })
  list(fits = fits, corrected = corr)
}
