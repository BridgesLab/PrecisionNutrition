# run_taskCD.R — Task C (sample overlap) + Task D (index-event / selection correction).
#
# Task C: quantify the exposure(CPD)-outcome(Bellenguez) UK Biobank overlap and document that
#   no UKB-excluded CPD GWAS is reachable in OpenGWAS (both ieu-b-142 and ieu-b-25 include UKB)
#   -> the overlap is instead broken from the OUTCOME side in Task B (Kunkle/Lambert have zero
#   UKB and the effect persists). MRlap requires full genome-wide sumstats + LD scores not
#   cached here -> reported as a reachable next step, not fabricated.
#
# Task D: SlopeHunter index-event correction. Selection (incidence) axis = parental lifespan;
#   prognosis = AD (Bellenguez). hunt() estimates the selection-induced slope b_SH. Because
#   IVW is linear, the corrected smoking->AD slope is exact:
#     corrected(smk->AD) = (smk->AD) - b_SH * (smk->lifespan)
#   applied to the pooled-nAChR and 15q25 estimates already computed (L6 + Task A).
#   Produces results/indexevent_correction.csv.

run_taskD <- function(cfg) {
  suppressPackageStartupMessages({library(ieugwasr); library(SlopeHunter); library(dplyr)})
  life_id <- cfg$selection_axis$parental_lifespan
  ad_id   <- cfg$opengwas$ad_primary

  # SNP set for the SlopeHunter fit: genome-wide tophits of the SELECTION (incidence) trait.
  th <- og_retry(function() tophits(life_id))            # clumped genome-wide-sig lifespan SNPs
  if (is.null(th) || nrow(th) < 10) return(list(status = "too few lifespan tophits", b_SH = NA))
  ad <- og_retry(function() extract_outcome_data(th$rsid, ad_id, proxies = TRUE, rsq = 0.8))
  if (is.null(ad) || nrow(ad) == 0) return(list(status = "no AD effects", b_SH = NA))

  dat <- th |>
    dplyr::transmute(SNP = rsid, BETA.incidence = beta, SE.incidence = se,
                     Pval.incidence = p, EA.incidence = ea, OA.incidence = nea) |>
    dplyr::inner_join(
      ad |> dplyr::transmute(SNP, BETA.prognosis = beta.outcome, SE.prognosis = se.outcome,
                             Pval.prognosis = pval.outcome,
                             EA.prognosis = effect_allele.outcome, OA.prognosis = other_allele.outcome),
      by = "SNP")
  # align prognosis effect alleles to incidence
  flip <- with(dat, EA.incidence != EA.prognosis & EA.incidence == OA.prognosis)
  dat$BETA.prognosis[flip] <- -dat$BETA.prognosis[flip]

  fit <- tryCatch(SlopeHunter::hunt(
    dat, snp_col = "SNP",
    xbeta_col = "BETA.incidence", xse_col = "SE.incidence", xp_col = "Pval.incidence",
    ybeta_col = "BETA.prognosis", yse_col = "SE.prognosis", yp_col = "Pval.prognosis",
    xp_thresh = 1, Plot = FALSE), error = function(e) {cat("hunt error:", conditionMessage(e), "\n"); NULL})
  if (is.null(fit)) return(list(status = "hunt failed", b_SH = NA, n_snp = nrow(dat)))
  b_SH <- as.numeric(fit$b)
  # SlopeHunter reports the slope SE under different field names across versions; fall back to
  # the bootstrap SD. With few selection-axis tophits this SE is large -> flagged downstream.
  se_SH <- suppressWarnings(as.numeric(fit$bse))
  if (length(se_SH) == 0 || is.na(se_SH)) se_SH <- stats::sd(fit$Bts.est, na.rm = TRUE)
  list(status = "ok", b_SH = b_SH, se_SH = se_SH, n_snp = nrow(dat),
       boot_sd = stats::sd(fit$Bts.est, na.rm = TRUE))
}

# Apply the linear index-event correction to a smoking->AD IVW estimate, given the matching
# smoking->lifespan IVW estimate and the SlopeHunter slope. Approx SE by delta method
# (treats the three estimates as independent — conservative-ish).
apply_correction <- function(b_ad, se_ad, b_life, se_life, b_SH, se_SH) {
  b_corr <- b_ad - b_SH * b_life
  se_corr <- sqrt(se_ad^2 + (b_SH^2) * se_life^2 + (b_life^2) * se_SH^2)
  tibble::tibble(b_uncorr = b_ad, se_uncorr = se_ad,
                 b_corr = b_corr, se_corr = se_corr,
                 ci_lo = b_corr - 1.96*se_corr, ci_hi = b_corr + 1.96*se_corr,
                 p_corr = 2*stats::pnorm(-abs(b_corr/se_corr)))
}
