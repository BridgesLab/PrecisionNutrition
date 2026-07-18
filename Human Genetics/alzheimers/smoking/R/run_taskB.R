# run_taskB.R — Task B: outcome ascertainment swap (the decisive test).
# Re-estimates the protective smoking->AD effect against a gradient of AD outcomes with
# DIFFERENT by-proxy fractions, on a harmonised per-SD-CPD log-OR scale, then tests whether
# the effect attenuates as proxy fraction falls (collider/selection signature) or persists
# (genuine effect). Sourced by 04c_ascertainment_swap.qmd; runnable standalone.
# Produces results/ascertainment_swap.csv and results/ascertainment_metareg.csv.

# Outcome gradient (composition-verified). proxy_frac = by-proxy cases / total cases.
# Bellenguez: 46828 proxy / (39106 clinical + 46828 proxy) = 0.545.
taskB_outcomes <- function(cfg) {
  a <- cfg$ad_ascertainment
  tibble::tribble(
    ~outcome,                ~id,              ~proxy_frac, ~scale,
    "Kunkle2019_clinical",   a$clinical_kunkle, 0.00,       "logOR",
    "Lambert2013_clinical",  a$clinical_lambert,0.00,       "logOR",
    "Bellenguez2022_proxymaj", a$proxy_majority, 0.545,     "logOR",
    "Schwartzentruber_GWAX", a$gwax_pure,       1.00,       "logOR_gwax")
}

# Standard (independent-SNP) IVW of the canonical 20-SNP set vs an outcome.
canonical20_mr <- function(outcome_id, cfg) {
  prn <- readr::read_csv(here::here("results","baseline_harmonised_pruned.csv"),
                         show_col_types = FALSE)
  exp <- prn |>
    dplyr::transmute(SNP, beta.exposure, se.exposure,
                     effect_allele.exposure, other_allele.exposure,
                     eaf.exposure, pval.exposure,
                     exposure = "smoking_CPD", id.exposure = "smoking_CPD",
                     mr_keep.exposure = TRUE, pval_origin.exposure = "reported")
  out <- og_retry(function() TwoSampleMR::extract_outcome_data(
    exp$SNP, outcome_id, proxies = TRUE, rsq = 0.8))
  if (is.null(out) || nrow(out) == 0) return(NULL)
  h <- suppressMessages(TwoSampleMR::harmonise_data(exp, out, action = 2)) |>
    dplyr::filter(mr_keep)
  if (nrow(h) < 2) return(NULL)
  m <- TwoSampleMR::mr(h, method_list = "mr_ivw")
  tibble::tibble(nsnp = m$nsnp[1], b = m$b[1], se = m$se[1], pval = m$pval[1])
}

# LD-aware cis IVW of a cached nAChR instrument set vs an outcome.
cis_set_mr <- function(instr_tag, outcome_id, label, cfg) {
  exp <- readRDS(here::here("data","cache",paste0(instr_tag,".rds")))
  mr <- ld_aware_cis_mr(exp, outcome_id, label, cfg)
  ivw <- mr |> dplyr::filter(method %in% c("LD-aware IVW","Wald ratio"))
  if (nrow(ivw) == 0 || is.na(ivw$b[1])) return(NULL)
  tibble::tibble(nsnp = ivw$nsnp[1], b = ivw$b[1], se = ivw$se[1], pval = ivw$pval[1])
}

# GWAX rescale: UKB family-history AD effects are ~half the case-control log-OR (each
# first-degree relative shares ~50% of genome). Multiply the GWAX MR estimate by 2 to
# approximate the case-control log-OR scale. Approximate; flagged in the writeup.
GWAX_SCALE <- 2

run_taskB <- function(cfg) {
  outs <- taskB_outcomes(cfg)
  sets <- tibble::tribble(
    ~set,             ~kind,    ~tag,
    "canonical_20snp","indep",  NA,
    "pooled_nAChR",   "cis",    "instr_pooled_nAChR",
    "15q25_CHRNA5A3B4","cis",   "instr_CHRNA5_A3_B4",
    "CHRNA4",         "cis",    "instr_CHRNA4")

  grid <- tidyr::expand_grid(set_i = seq_len(nrow(sets)), out_i = seq_len(nrow(outs)))
  res <- purrr::pmap_dfr(grid, function(set_i, out_i) {
    s <- sets[set_i,]; o <- outs[out_i,]
    est <- if (s$kind == "indep") canonical20_mr(o$id, cfg)
           else cis_set_mr(s$tag, o$id, s$set, cfg)
    if (is.null(est))
      return(tibble::tibble(set = s$set, outcome = o$outcome, id = o$id,
                            proxy_frac = o$proxy_frac, nsnp = 0L,
                            b = NA_real_, se = NA_real_, pval = NA_real_, status = "MR failed"))
    sc <- if (o$scale == "logOR_gwax") GWAX_SCALE else 1
    tibble::tibble(set = s$set, outcome = o$outcome, id = o$id,
                   proxy_frac = o$proxy_frac, nsnp = est$nsnp,
                   b = est$b * sc, se = est$se * sc, pval = est$pval,
                   scale = o$scale, status = "ok")
  })
  res |> dplyr::mutate(or = exp(b), ci_lo = b - 1.96*se, ci_hi = b + 1.96*se)
}

# Proxy->clinical difference test + meta-regression of effect on proxy fraction, per set.
# Clinical anchor = inverse-variance mean of the two clinical (proxy=0) outcomes.
taskB_tests <- function(res) {
  ivw_pool <- function(b, se) {
    w <- 1/se^2; bm <- sum(w*b)/sum(w); sem <- sqrt(1/sum(w))
    c(b = bm, se = sem)
  }
  out <- res |> dplyr::filter(status == "ok") |> dplyr::group_split(set) |>
    purrr::map_dfr(function(d) {
      set <- d$set[1]
      clin <- d |> dplyr::filter(proxy_frac == 0)
      prox <- d |> dplyr::filter(outcome == "Bellenguez2022_proxymaj")
      zdiff <- NA; pdiff <- NA; b_clin <- NA; se_clin <- NA
      if (nrow(clin) >= 1 && nrow(prox) == 1) {
        pc <- ivw_pool(clin$b, clin$se); b_clin <- pc["b"]; se_clin <- pc["se"]
        zdiff <- (prox$b - b_clin) / sqrt(prox$se^2 + se_clin^2)
        pdiff <- 2*pnorm(-abs(zdiff))
      }
      # meta-regression b ~ proxy_frac, inverse-variance weighted, across all outcomes
      slope <- NA; slope_se <- NA; slope_p <- NA
      if (nrow(d) >= 3) {
        fit <- tryCatch(stats::lm(b ~ proxy_frac, data = d, weights = 1/se^2),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          co <- summary(fit)$coefficients
          if ("proxy_frac" %in% rownames(co)) {
            slope <- co["proxy_frac","Estimate"]; slope_se <- co["proxy_frac","Std. Error"]
            slope_p <- co["proxy_frac","Pr(>|t|)"]
          }
        }
      }
      tibble::tibble(set = set, b_clinical = unname(b_clin), se_clinical = unname(se_clin),
                     b_proxymaj = if (nrow(prox)==1) prox$b else NA,
                     z_diff = zdiff, p_diff = pdiff,
                     metareg_slope = slope, metareg_se = slope_se, metareg_p = slope_p)
    })
  out
}
