# mvmr_lifespan.R — Zhang et al. 2026's step 2 for the survival collider: multivariable MR with the
# exposure and lifespan as joint exposures on the RAW AD outcome.
#
#   beta_AD(j) = theta_x * beta_x(j) + theta_l * beta_life(j) + error
#
# The AD GWAS is conditioned on survival, so each variant's AD association carries b * beta_life
# (Dudbridge 2019; Zhang et al. 2026, eq. 2). theta_l absorbs that slope (compare b_SH = -0.944) and
# theta_x is the exposure effect with it removed: the same correction as the SlopeHunter arm, with b
# estimated from lifespan instruments instead of SlopeHunter's mixture model, and its uncertainty
# carried by the regression. Nothing is conditioned on individually; lifespan enters only through
# genetic associations estimated in the (unselected) lifespan GWAS.
#
# Validity needs lifespan instruments that affect AD only through survival. Zhang et al. (Supp.
# Table S6) show the bias when they do not, so lifespan instruments at established AD loci are
# excluded (MVMR_EXCLUDE); the sensitivity set keeps all but APOE.
#
# Sources pipeline/R helpers: plink_clump(), error_cor(), wls_mre(), mrbee_mv(), fit_rows(), try_fit().

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# hg19 windows. APOE matches the SlopeHunter fits (it spans both builds' APOE coordinates).
MVMR_EXCLUDE <- tibble::tribble(
  ~region,     ~chr, ~start, ~end,   ~why,
  "APOE",      "19", 43.9e6, 46.5e6, "AD locus; excluded from every SlopeHunter fit too",
  "PTK2B-CLU", "8",  27.1e6, 27.7e6, "established AD loci (Lambert 2013; Bellenguez 2022)",
  "MHC",       "6",  28.4e6, 33.5e6, "HLA-DQA1/DRB1, established AD locus (Lambert 2013; Bellenguez 2022)")

in_region <- function(chr, pos, regions = MVMR_EXCLUDE) {
  hit <- rep(NA_character_, length(chr))
  for (k in seq_len(nrow(regions))) {
    w <- as.character(chr) == regions$chr[k] & pos > regions$start[k] & pos < regions$end[k]
    hit[w & is.na(hit)] <- regions$region[k]
  }
  hit
}


# Genome-wide significant, independent lifespan variants, flagged by exclusion region.
life_instruments <- function(g_life, bfile, plink2, p = 5e-8, r2 = 0.001, kb = 10000) {
  idx <- plink_clump(g_life, bfile, p, r2, kb, plink2)
  iv <- g_life[SNP %chin% idx]
  iv[, region := in_region(chr, pos)]
  iv[order(p)]
}

# g's beta/se/p on ref's effect allele, through the pipeline's harmonise_pair() (TwoSampleMR,
# action = 2: palindromes resolved by allele frequency, ambiguous ones dropped). Columns suffixed.
align_to <- function(ref, g, sfx) {
  inst <- as.data.table(ref)[, .(SNP, chr, pos, ea, oa, eaf, beta, se, p, n)]
  h <- harmonise_pair(inst, g, "mvmr_lifespan", "mvmr")$dat
  out <- if (is.null(h)) data.table(SNP = character(), beta = numeric(), se = numeric(), p = numeric())
         else as.data.table(h)[mr_keep %in% TRUE, .(SNP, beta = beta.outcome, se = se.outcome, p = pval.outcome)]
  setnames(out, c("beta", "se", "p"), paste0(c("beta_", "se_", "p_"), sfx))
  out
}

# Joint instrument set: exposure instruments + the lifespan instruments in use, APOE region removed,
# re-clumped with LIFESPAN PRIORITY. Clumping on the smaller p-value (the usual construction) lets
# hundreds of strong body-composition instruments win nearly every LD conflict (r2 0.001 over 10 Mb),
# leaving 1-3 lifespan variants and a lifespan conditional F of ~2, so the correction is not
# identified. Same fix as mr_mediation.qmd's mediator-priority clumping. Lifespan variants missing
# from the exposure GWAS (the HapMap-era fasting-glucose files) are lost; no proxies are sought.
# Rows are the exposure GWAS's (its alleles are the reference).
mvmr_instruments <- function(inst_x, g_x, g_life, life_snps, bfile, plink2, r2 = 0.001, kb = 10000) {
  ids <- union(inst_x$SNP, life_snps)
  u <- merge(g_x[SNP %chin% ids], g_life[SNP %chin% ids, .(SNP, p_l = p)], by = "SNP", all.x = TRUE)
  u[, `:=`(p_x = p, source = fifelse(SNP %chin% life_snps, "lifespan", "exposure"))]
  u[, p := fifelse(source == "lifespan", 1e-300, p_x)]
  u <- u[!in_region(chr, pos, MVMR_EXCLUDE[1, ]) %in% "APOE"]
  idx <- plink_clump(u, bfile, 5e-8, r2, kb, plink2)
  u[SNP %chin% idx][, p := p_x][]
}

# Exposure, lifespan and outcome effects on the exposure's effect allele, for the variants all three
# GWAS can harmonise.
mvmr_data <- function(ref, g_life, g_y)
  Reduce(\(a, b) merge(a, b, by = "SNP"),
         list(as.data.table(ref)[, .(SNP, chr, pos, source, beta_x = beta, se_x = se, p_x)],
              align_to(ref, g_life, "l"), align_to(ref, g_y, "y")))

# Error-correlation matrix (exposure, lifespan, outcome) for MRBEE.
mvmr_R <- function(r) matrix(c(1, r[["x_l"]], r[["x_y"]],
                               r[["x_l"]], 1, r[["l_y"]],
                               r[["x_y"]], r[["l_y"]], 1), 3)

# MVMR by IVW-MRE and MRBEE, plus the univariable fit on the exposure-selected variants of the same
# data as the uncorrected comparator.
mvmr_fit <- function(d, r) {
  X <- cbind(exposure = d$beta_x, lifespan = d$beta_l)
  S <- cbind(d$se_x, d$se_l)
  u <- d[d$source == "exposure", ]
  Xu <- cbind(exposure = u$beta_x)
  Ru <- matrix(c(1, r[["x_y"]], r[["x_y"]], 1), 2)
  bind_rows(
    try_fit(fit_rows(wls_mre(d$beta_y, X, d$se_y), n = nrow(d)), "MVMR", "IVW-MRE"),
    try_fit(fit_rows(mrbee_mv(d$beta_y, X, d$se_y, S, mvmr_R(r)), n = nrow(d)), "MVMR", "MRBEE"),
    try_fit(fit_rows(wls_mre(u$beta_y, Xu, u$se_y), n = nrow(u)), "univariable", "IVW-MRE"),
    try_fit(fit_rows(mrbee_uv(u$beta_y, u$beta_x, u$se_y, u$se_x, Ru), n = nrow(u)), "univariable", "MRBEE"))
}

# Univariable MRBEE in fit_rows() form (as classical.R calls it).
mrbee_uv <- function(y, x, se_y, se_x, R) {
  f <- suppressWarnings(MRBEE::MRBEE.IMRP.UV(by = y, bx = x, byse = se_y, bxse = se_x, Rxy = R))
  list(b = c(exposure = as.numeric(f$theta)[1]),
       V = matrix(as.numeric(f$vartheta)[1], 1, 1, dimnames = list("exposure", "exposure")),
       sigma = NA_real_, n_outlier = if (!is.null(f$delta)) sum(abs(as.numeric(f$delta)) > 1e-8) else NA_integer_)
}

# Sanderson-Windmeijer conditional F for both exposures (MVMR::strength_mvmr), with the
# exposure-lifespan error covariance from sample overlap.
mvmr_cond_f <- function(d, r_xl) {
  se <- cbind(d$se_x, d$se_l)
  inp <- MVMR::format_mvmr(BXGs = cbind(d$beta_x, d$beta_l), BYG = d$beta_y, seBXGs = se,
                           seBYG = d$se_y, RSID = d$SNP)
  gencov <- MVMR::phenocov_mvmr(matrix(c(1, r_xl, r_xl, 1), 2), se)
  utils::capture.output(s <- MVMR::strength_mvmr(inp, gencov))
  c(exposure = s[1, 1], lifespan = s[1, 2])
}

# The selection slope from lifespan instruments alone (AD on lifespan): IVW-MRE, weighted median,
# MR-Egger. This is what theta_l estimates, without an exposure in the model.
life_slope <- function(iv, g_y) {
  d <- merge(iv[, .(SNP, ea, oa, beta_l = beta, se_l = se)], align_to(iv, g_y, "y"), by = "SNP")
  f <- wls_mre(d$beta_y, cbind(lifespan = d$beta_l), d$se_y)
  wm <- TwoSampleMR::mr_weighted_median(d$beta_l, d$beta_y, d$se_l, d$se_y)
  eg <- TwoSampleMR::mr_egger_regression(d$beta_l, d$beta_y, d$se_l, d$se_y)
  tibble(method = c("IVW-MRE", "Weighted median", "MR-Egger"),
         b = c(f$b[[1]], wm$b, eg$b), se = c(sqrt(f$V[1, 1]), wm$se, eg$se),
         egger_intercept_p = c(NA, NA, eg$pval_i), n_iv = nrow(d),
         Q_p = stats::pchisq(f$sigma^2 * (nrow(d) - 1), nrow(d) - 1, lower.tail = FALSE)) |>
    mutate(p = 2 * stats::pnorm(-abs(b / se)), .after = se)
}

# Leave one lifespan-selected variant out at a time (IVW-MRE and MRBEE).
mvmr_loo <- function(d, r) {
  life <- d$SNP[d$source == "lifespan"]
  purrr::map_dfr(life, \(s) mvmr_fit(d[d$SNP != s, ], r) |> filter(step == "MVMR") |>
                   mutate(dropped = s, .before = 1))
}

# Driver: every exposure x outcome x lifespan-instrument set. GWAS come from the targets store
# (gwas_<key>, inst_classical_<key>). Writes CSVs to out_dir and returns them as a list.
run_mvmr_lifespan <- function(exposures, outcomes, bfile, plink2, out_dir = "results/mvmr_lifespan") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  rd <- \(k) targets::tar_read_raw(k)
  g_life <- rd("gwas_lifespan")
  g_y <- purrr::set_names(purrr::map(outcomes, \(o) rd(paste0("gwas_", o))), outcomes)
  iv <- life_instruments(g_life, bfile, plink2)
  sets <- list(primary = iv$SNP[is.na(iv$region)], all_but_apoe = iv$SNP[!iv$region %in% "APOE"])
  r_ly <- purrr::map_dbl(g_y, \(g) error_cor(g_life, g)[1, 2])

  iv_out <- purrr::imap_dfr(g_y, \(g, o) merge(iv, align_to(iv, g, "y"), by = "SNP", all.x = TRUE) |>
                              as_tibble() |> mutate(outcome = o)) |>
    select(outcome, SNP, chr, pos, ea, oa, region, beta_l = beta, se_l = se, p_l = p, beta_y, se_y, p_y)
  slope <- purrr::imap_dfr(sets, \(snps, set) purrr::imap_dfr(g_y, \(g, o)
    life_slope(iv[SNP %chin% snps], g) |> mutate(outcome = o, set = set, .before = 1)))

  res <- list(); cf <- list(); loo <- list(); dat <- list()
  for (e in exposures) {
    message("MVMR: ", e)
    g_x <- rd(paste0("gwas_", e)); inst_x <- rd(paste0("inst_classical_", e))$inst
    r_xl <- error_cor(g_x, g_life)[1, 2]
    r_xy <- purrr::map_dbl(g_y, \(g) error_cor(g_x, g)[1, 2])
    for (set in names(sets)) {
      ref <- mvmr_instruments(inst_x, g_x, g_life, sets[[set]], bfile, plink2)
      for (o in outcomes) {
        d <- as_tibble(mvmr_data(ref, g_life, g_y[[o]]))
        r <- c(x_l = r_xl, x_y = r_xy[[o]], l_y = r_ly[[o]])
        tag <- \(x) mutate(x, exposure = e, outcome = o, set = set, .before = 1)
        res[[length(res) + 1]] <- tag(mvmr_fit(d, r)) |>
          mutate(n_life_snp = sum(d$source == "lifespan"), r_xl = r_xl, r_xy = r[["x_y"]], r_ly = r[["l_y"]])
        f <- tryCatch(mvmr_cond_f(d, r_xl), error = \(err) c(exposure = NA_real_, lifespan = NA_real_))
        cf[[length(cf) + 1]] <- tag(tibble(n_snp = nrow(d), n_life_snp = sum(d$source == "lifespan"),
                                           cond_F_exposure = f[["exposure"]], cond_F_lifespan = f[["lifespan"]]))
        if (set == "primary") {
          loo[[length(loo) + 1]] <- tag(mvmr_loo(d, r))
          dat[[length(dat) + 1]] <- tag(d)
        }
      }
    }
    rm(g_x); gc(verbose = FALSE)
  }
  out <- list(life_instruments = iv_out, life_slope = slope, results = bind_rows(res),
              cond_f = bind_rows(cf), loo = bind_rows(loo), data = bind_rows(dat))
  purrr::iwalk(out, \(x, nm) readr::write_csv(x, file.path(out_dir, paste0(nm, ".csv"))))
  out
}
