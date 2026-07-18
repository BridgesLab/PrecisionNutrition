# run_taskA.R — Task A: locus-resolved selection negative-control battery.
# Tests the EFFECT-DRIVING nAChR instrument sets (not just CHRNB2) against the
# selection/survival axis with LD-aware IVW + pre-registered TOST equivalence (APPROACH.md
# sec 13). Sourced by 04b_selection_robustness.qmd and runnable standalone via Rscript.
# Produces results/selection_battery.csv. Instruments are the cached L6 cis sets (fixed).

run_taskA <- function(cfg, delta_primary = 0.02, delta_sens = 0.05) {
  # Instrument sets (cached by get_locus_instruments / pooled in 04b setup).
  instr_sets <- list(
    `15q25_CHRNA5A3B4` = "instr_CHRNA5_A3_B4",
    CHRNA4             = "instr_CHRNA4",
    CHRNB2             = "instr_CHRNB2",
    CHRNA6B3           = "instr_CHRNA6_B3",
    pooled_nAChR       = "instr_pooled_nAChR")

  # Selection-axis outcomes (composition-verified; see dataset_composition_verified.csv).
  sa <- cfg$selection_axis
  outcomes <- tibble::tribble(
    ~outcome,            ~id,                ~scale,        ~ukb,
    "parental_lifespan", sa$parental_lifespan, "SD/quant",  TRUE,
    "lung_cancer",       sa$lung_cancer,       "log-OR",    FALSE,
    "COPD_FEV1FVC",      sa$copd_fev1fvc,      "SD/quant",  TRUE,
    "CAD",               sa$cad_primary,       "log-OR",    FALSE)

  load_instr <- function(tag) {
    f <- here::here("data", "cache", paste0(tag, ".rds"))
    if (file.exists(f)) readRDS(f) else NULL
  }

  grid <- tidyr::expand_grid(set = names(instr_sets), outcome = outcomes$outcome)
  res <- purrr::pmap_dfr(grid, function(set, outcome) {
    exp <- load_instr(instr_sets[[set]])
    oid <- outcomes$id[outcomes$outcome == outcome]
    if (is.null(exp) || is.na(oid)) {
      return(tibble::tibble(set = set, outcome = outcome, id = oid %||% NA, nsnp = 0L,
                            b = NA_real_, se = NA_real_, status = "no instrument/outcome"))
    }
    mr <- ld_aware_cis_mr(exp, oid, set, cfg)
    ivw <- mr |> dplyr::filter(method %in% c("LD-aware IVW", "Wald ratio"))
    if (nrow(ivw) == 0 || is.na(ivw$b[1]))
      return(tibble::tibble(set = set, outcome = outcome, id = oid, nsnp = 0L,
                            b = NA_real_, se = NA_real_, status = "MR failed"))
    tibble::tibble(set = set, outcome = outcome, id = oid, nsnp = ivw$nsnp[1],
                   method = ivw$method[1], b = ivw$b[1], se = ivw$se[1], status = "ok")
  })

  # Attach TOST equivalence at both bounds + the pre-registered verdict.
  eq_p <- purrr::pmap_dfr(list(res$b, res$se), function(b, s) tost_equivalence(b, s, delta_primary))
  eq_s <- purrr::pmap_dfr(list(res$b, res$se), function(b, s) tost_equivalence(b, s, delta_sens))
  res |>
    dplyr::mutate(or = exp(b),
                  p_assoc = eq_p$p_assoc, ci_lo = eq_p$ci_lo, ci_hi = eq_p$ci_hi,
                  tost_p_d02 = eq_p$tost_p, verdict_d02 = eq_p$verdict,
                  tost_p_d05 = eq_s$tost_p, verdict_d05 = eq_s$verdict) |>
    dplyr::left_join(dplyr::select(outcomes, outcome, scale, ukb), by = "outcome")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Task A, Part 2: the Y->S arm of the survival collider ---------------------------------
# Does AD liability feed the survival/selection node? Estimates AD -> parental lifespan by MR.
# Uses CLINICAL AD instruments (Kunkle ieu-b-2) as the exposure -- NEVER a proxy AD GWAS, which
# would be circular (proxy-AD is itself built on parental survival). APOE dominates AD->lifespan
# via genuine pleiotropy (ApoE4 -> AD AND -> shorter life) and is NOT in the nAChR instrument
# set, so we report the estimate with and without the APOE region: the non-APOE AD axis is the
# part relevant to whether the collider operates at the receptor loci.
# Produces results/collider_Yarm_AD_to_lifespan.csv.

# Well-known APOE-region tag SNPs (GRCh37 chr19:~44.4-45.5Mb), used as a fallback when the
# harmonised frame lacks position columns.
APOE_SNPS <- c("rs429358","rs7412","rs4420638","rs769449","rs6857","rs2075650",
               "rs439401","rs445925","rs56131196","rs157582","rs59007384")

run_taskA_Yarm <- function(cfg) {
  # Clinical AD instrument sources, tried in order (never a proxy GWAS — circular).
  ad_sources <- c(cfg$ad_ascertainment$clinical_kunkle,    # ieu-b-2 Kunkle 2019
                  cfg$ad_ascertainment$clinical_lambert)   # ieu-a-297 Lambert 2013 (fallback)
  life_id <- cfg$selection_axis$parental_lifespan          # Pilling combined parental attained age

  # Consistent empty row (full schema) so the qmd table renders even with no reachable data.
  empty_y <- function(status) tibble::tibble(
    snp_set = status, nsnp = 0L, b_ivw = NA_real_, se_ivw = NA_real_, p_ivw = NA_real_,
    b_wmed = NA_real_, se_wmed = NA_real_, p_wmed = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, status = status)

  inst <- NULL; ad_id <- NA
  for (id in ad_sources) {
    cand <- og_retry(function() TwoSampleMR::extract_instruments(id, p1 = 5e-8, clump = TRUE))
    if (!is.null(cand) && nrow(cand) > 0) { inst <- cand; ad_id <- id; break }
  }
  if (is.null(inst) || nrow(inst) == 0) return(empty_y("no clinical AD instrument reachable"))
  message("Y-arm AD instrument source: ", ad_id)
  out <- og_retry(function() TwoSampleMR::extract_outcome_data(inst$SNP, life_id,
                                                               proxies = TRUE, rsq = 0.8))
  if (is.null(out) || nrow(out) == 0) return(empty_y("no lifespan outcome"))
  h <- suppressMessages(TwoSampleMR::harmonise_data(inst, out, action = 2)) |>
    dplyr::filter(mr_keep)
  if (nrow(h) < 2) return(empty_y("too few harmonised SNP"))

  # APOE-region flag: by position if available, else by tag-SNP fallback.
  pos_col <- intersect(c("pos.exposure","position","pos"), names(h))
  chr_col <- intersect(c("chr.exposure","chr"), names(h))
  if (length(pos_col) && length(chr_col)) {
    apoe <- as.integer(h[[chr_col[1]]]) == 19 &
            as.numeric(h[[pos_col[1]]]) > 44.0e6 & as.numeric(h[[pos_col[1]]]) < 46.5e6
    apoe[is.na(apoe)] <- h$SNP[is.na(apoe)] %in% APOE_SNPS
  } else {
    apoe <- h$SNP %in% APOE_SNPS
  }

  est_one <- function(hh, label) {
    if (nrow(hh) < 2) return(NULL)
    m <- TwoSampleMR::mr(hh, method_list = c("mr_ivw","mr_weighted_median"))
    ivw <- m |> dplyr::filter(method == "Inverse variance weighted")
    wm  <- m |> dplyr::filter(method == "Weighted median")
    tibble::tibble(snp_set = label, nsnp = ivw$nsnp[1],
                   b_ivw = ivw$b[1], se_ivw = ivw$se[1], p_ivw = ivw$pval[1],
                   b_wmed = wm$b[1], se_wmed = wm$se[1], p_wmed = wm$pval[1],
                   ci_lo = ivw$b[1] - 1.96*ivw$se[1], ci_hi = ivw$b[1] + 1.96*ivw$se[1],
                   status = "ok")
  }
  dplyr::bind_rows(
    est_one(h, sprintf("all AD instruments (n=%d)", nrow(h))),
    est_one(h[!apoe, ], sprintf("excl APOE region (n=%d)", sum(!apoe)))
  )
}
