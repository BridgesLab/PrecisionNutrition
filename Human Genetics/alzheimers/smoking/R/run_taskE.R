# run_taskE.R — Task E: direction negative controls on the OUTCOME itself.
# Independent of any nAChR claim: does the primary AD outcome (Bellenguez, proxy-majority)
# behave like a survival/participation-selected phenotype? A selected outcome makes
# education appear paradoxically protective and lung cancer appear "protective" against AD.
# Runs EA->AD and lung-cancer->AD (and CAD->AD) through the SAME primary outcome.
# Produces results/direction_negcontrols.csv.

run_taskE <- function(cfg) {
  exposures <- tibble::tribble(
    ~exposure,       ~id,                ~note,
    "education",     "ieu-a-1239",        "Lee 2018 EA (excl 23andMe)",
    "lung_cancer",   "ebi-a-GCST004748",  "McKay 2017 ILCCO",
    "CAD",           "ieu-a-7",           "Nikpay 2015")
  ad <- cfg$opengwas$ad_primary

  one <- function(exposure, id, note) {
    inst <- og_retry(function() TwoSampleMR::extract_instruments(id, p1 = 5e-8, clump = TRUE))
    if (is.null(inst) || nrow(inst) == 0)
      return(tibble::tibble(exposure = exposure, id = id, nsnp = 0L,
                            b = NA_real_, se = NA_real_, pval = NA_real_, status = "no instrument"))
    out <- og_retry(function() TwoSampleMR::extract_outcome_data(inst$SNP, ad, proxies = TRUE, rsq = 0.8))
    if (is.null(out) || nrow(out) == 0)
      return(tibble::tibble(exposure = exposure, id = id, nsnp = 0L,
                            b = NA_real_, se = NA_real_, pval = NA_real_, status = "no outcome"))
    h <- suppressMessages(TwoSampleMR::harmonise_data(inst, out, action = 2)) |>
      dplyr::filter(mr_keep)
    if (nrow(h) < 2)
      return(tibble::tibble(exposure = exposure, id = id, nsnp = nrow(h),
                            b = NA_real_, se = NA_real_, pval = NA_real_, status = "too few SNP"))
    m <- TwoSampleMR::mr(h, method_list = c("mr_ivw","mr_weighted_median"))
    ivw <- m |> dplyr::filter(method == "Inverse variance weighted")
    tibble::tibble(exposure = exposure, id = id, note = note, nsnp = ivw$nsnp[1],
                   b = ivw$b[1], se = ivw$se[1], pval = ivw$pval[1],
                   or = exp(ivw$b[1]), ci_lo = ivw$b[1]-1.96*ivw$se[1], ci_hi = ivw$b[1]+1.96*ivw$se[1],
                   status = "ok")
  }
  purrr::pmap_dfr(exposures, one) |>
    dplyr::mutate(
      # Selection signature: lung cancer / CAD "protective" against AD, education paradoxical.
      selection_signature = dplyr::case_when(
        exposure %in% c("lung_cancer","CAD") & !is.na(b) & b < 0 & pval < 0.05 ~ "protective (SELECTION SIGNATURE)",
        exposure == "education" & !is.na(b) & b < 0 & pval < 0.05 ~ "protective (expected if outcome selected/proxy)",
        !is.na(b) & pval >= 0.05 ~ "null",
        TRUE ~ "other"))
}
