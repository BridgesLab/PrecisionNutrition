# classical.R — the full conventional MR suite and instrument diagnostics for one association.
#
# Runs on the classical instrument set (p < 5e-8, r2 < 0.001). Estimates are per unit of each GWAS's
# own beta (log-OR for binary traits), unlike MRAID/CAUSE, which work on standardised effects.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# Correlation of exposure/outcome z-scores among variants null for both (p > p_min): the
# estimation-error correlation MRBEE needs. Non-zero mostly through sample overlap. Same estimator
# as calcium-cholesterol/R/robust_mr_helpers.R::mrbee_error_cov(), computed genome-wide here.
error_cor <- function(g_exp, g_out, p_min = 0.05, max_n = 500000, seed = 1) {
  m <- merge(g_exp[, .(SNP, ea, oa, zx = beta / se, px = p)],
             g_out[, .(SNP, ea_o = ea, oa_o = oa, zy = beta / se, py = p)], by = "SNP")
  m <- m[(ea == ea_o & oa == oa_o) | (ea == oa_o & oa == ea_o)]
  m[ea != ea_o, zy := -zy]
  m <- m[px > p_min & py > p_min]
  if (nrow(m) > max_n) { set.seed(seed); m <- m[sample(.N, max_n)] }
  R <- stats::cor(cbind(m$zx, m$zy))
  dimnames(R) <- list(c("exposure", "outcome"), c("exposure", "outcome"))
  R
}

# Instrument-strength and heterogeneity diagnostics.
#   Q, I2      Cochran's Q for IVW and I2 = (Q - df)/Q
#   Q_rucker   Rucker's Q' from MR-Egger (Q - Q' > 0 favours Egger's intercept)
#   I2_GX      Bowden 2016: precision of the SNP-exposure effects; < 0.9 means Egger is diluted
#   R2         sum over SNPs of z^2 / (z^2 + N - 2), z = beta_x/se_x (scale-free, so it is valid for
#              the raw-mmol/L 2hGlu betas; observed scale for binary exposures)
#   F          mean of per-SNP z^2; F_total = R2 (N - k - 1) / ((1 - R2) k)
instrument_diagnostics <- function(dat) {
  bx <- dat$beta.exposure; sx <- dat$se.exposure; by <- dat$beta.outcome; sy <- dat$se.outcome
  k <- nrow(dat); n <- stats::median(dat$samplesize.exposure, na.rm = TRUE)
  w <- 1 / sy^2
  b_ivw <- sum(w * bx * by) / sum(w * bx^2)
  Q <- sum(w * (by - b_ivw * bx)^2)
  het <- tryCatch(TwoSampleMR::mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw")),
                  error = function(e) NULL)
  q_rucker <- if (!is.null(het)) het$Q[het$method == "MR Egger"] else NA_real_
  y <- abs(bx); ws <- 1 / sx^2; ybar <- sum(ws * y) / sum(ws)
  Qgx <- sum(ws * (y - ybar)^2)
  z2 <- (bx / sx)^2
  r2 <- if (is.finite(n)) sum(z2 / (z2 + n - 2)) else NA_real_
  tibble(n_snp = k, Q = Q, Q_df = k - 1, Q_p = stats::pchisq(Q, k - 1, lower.tail = FALSE),
         I2 = max(0, (Q - (k - 1)) / Q),
         Q_rucker = q_rucker %||% NA_real_,
         Q_rucker_p = if (!is.null(het)) het$Q_pval[het$method == "MR Egger"] else NA_real_,
         I2_GX = max(0, (Qgx - (k - 1)) / Qgx),
         R2 = r2, N_exposure = n, mean_F = mean(z2),
         total_F = if (is.finite(r2)) r2 * (n - k - 1) / ((1 - r2) * k) else NA_real_)
}

# Steiger directionality with the same z-based R2 on both sides (approximate for binary traits, which
# have no liability-scale conversion here; flagged in the output).
steiger_z <- function(dat) {
  r2 <- function(b, s, n) { z2 <- (b / s)^2; z2 / (z2 + n - 2) }
  rx <- r2(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
  ry <- r2(dat$beta.outcome, dat$se.outcome, dat$samplesize.outcome)
  tibble(steiger_r2_exposure = sum(rx, na.rm = TRUE), steiger_r2_outcome = sum(ry, na.rm = TRUE),
         steiger_correct_direction = sum(rx, na.rm = TRUE) > sum(ry, na.rm = TRUE),
         steiger_n_snp_wrong = sum(ry > rx, na.rm = TRUE))
}

# One method -> one row. Failures are rows with an error message, never silent drops.
method_row <- function(method, b = NA_real_, se = NA_real_, p = NA_real_, n = NA_integer_, ...)
  tibble(method = method, b = as.numeric(b), se = as.numeric(se),
         p = if (is.na(p) && is.finite(b) && is.finite(se) && se > 0) 2 * pnorm(-abs(b / se)) else as.numeric(p),
         n_snp = as.integer(n), ...)

run_classical_suite <- function(dat, Rxy = NULL, presso_nboot = 1000, seed = 20261004) {
  if (is.null(dat) || nrow(dat) < 3)
    return(method_row("all", note = sprintf("not run: %d instrument(s)", nrow(dat %||% data.frame()))))
  set.seed(seed)
  tsm <- function(code, label) {
    r <- tryCatch(suppressMessages(TwoSampleMR::mr(dat, method_list = code)), error = function(e) e)
    if (inherits(r, "error") || !nrow(r)) return(method_row(label, note = paste("failed:", conditionMessage(r %||% simpleError("no result")))))
    method_row(label, r$b[1], r$se[1], r$pval[1], r$nsnp[1])
  }
  # IVW-MRE computed here, not with TwoSampleMR::mr_ivw_mre(): in TwoSampleMR 0.7.5 that function
  # uses the raw residual SE with no floor, so under-dispersed instruments (Q < df) get an SE
  # SMALLER than fixed effects. The standard multiplicative random-effects model (Bowden et al.
  # 2017; MendelianRandomization::mr_ivw(model = "random")) scales by max(1, residual SE).
  ivw_mre <- function() {
    s <- summary(stats::lm(dat$beta.outcome ~ -1 + dat$beta.exposure, weights = 1 / dat$se.outcome^2))
    b <- s$coefficients[1, 1]
    method_row("IVW-MRE", b, s$coefficients[1, 2] / s$sigma * max(1, s$sigma), NA, nrow(dat),
               ivw_residual_se = s$sigma)
  }
  rows <- list(
    ivw_mre(), tsm("mr_ivw_fe", "IVW-FE"), tsm("mr_egger_regression", "MR-Egger"),
    tsm("mr_weighted_median", "Weighted median"), tsm("mr_weighted_mode", "Weighted mode"),
    tsm("mr_raps", "MR-RAPS"))

  egg <- tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error = function(e) NULL)
  rows[[3]] <- rows[[3]] |> mutate(egger_intercept = egg$egger_intercept %||% NA_real_,
                                   egger_intercept_se = egg$se %||% NA_real_,
                                   egger_intercept_p = egg$pval %||% NA_real_)

  pr <- if (nrow(dat) >= 4) tryCatch(suppressWarnings(MRPRESSO::mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
    SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = as.data.frame(dat), NbDistribution = presso_nboot, SignifThreshold = 0.05)),
    error = function(e) e) else simpleError("fewer than 4 SNPs")
  rows[[length(rows) + 1]] <- if (inherits(pr, "error")) method_row("MR-PRESSO", note = paste("failed:", conditionMessage(pr))) else {
    num <- function(x) suppressWarnings(as.numeric(sub("^<", "", as.character(x))))
    mm <- pr$`Main MR results`
    oc <- mm[mm$`MR Analysis` == "Outlier-corrected", ]
    use <- if (nrow(oc) && is.finite(num(oc$`Causal Estimate`))) oc else mm[mm$`MR Analysis` == "Raw", ]
    idx <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    method_row("MR-PRESSO", num(use$`Causal Estimate`), num(use$Sd), num(use$`P-value`),
               nrow(dat) - if (is.numeric(idx)) length(idx) else 0L,
               presso_global_p = num(pr$`MR-PRESSO results`$`Global Test`$Pvalue),
               presso_n_outlier = if (is.numeric(idx)) length(idx) else 0L,
               presso_estimate = use$`MR Analysis`[1],
               presso_distortion_p = num(pr$`MR-PRESSO results`$`Distortion Test`$Pvalue %||% NA))
  }

  mb <- if (is.null(Rxy)) simpleError("no error-correlation matrix") else tryCatch(suppressWarnings(
    MRBEE::MRBEE.IMRP.UV(by = dat$beta.outcome, bx = dat$beta.exposure, byse = dat$se.outcome,
                         bxse = dat$se.exposure, Rxy = Rxy)), error = function(e) e)
  rows[[length(rows) + 1]] <- if (inherits(mb, "error")) method_row("MRBEE", note = paste("failed:", conditionMessage(mb))) else
    method_row("MRBEE", as.numeric(mb$theta)[1], sqrt(as.numeric(mb$vartheta)[1]), NA, nrow(dat),
               mrbee_n_pleiotropic = if (!is.null(mb$delta)) sum(abs(as.numeric(mb$delta)) > 1e-8) else NA_integer_,
               mrbee_error_cor = Rxy[1, 2])
  bind_rows(rows) |> mutate(scale = "per unit of GWAS beta", .after = method)
}
