# contrast_outcome.R — fasting vs 2 h glucose -> AD, Steps 1-3 and the drug-target layer.
#
# UNBLINDED: everything here reads an AD GWAS. Each entry point calls assert_frozen() first, so it
# cannot run unless assignments/ matches its manifest and is committed (the Step 0 freeze).
#
# All SNP effects are on the T2D risk allele. AD effects are log-odds. Two exposure scales:
#   mmol   betas as in the MAGIC files (mmol/L)                       -> log-OR AD per mmol/L
#   sd     betas divided by the trait SD (config contrast: sd_mmol)   -> log-OR AD per SD
#
# Step 1  MVMR        beta_AD ~ th_fg * b_fg + th_2h * b_2h; headline = Wald test of th_2h - th_fg
#                     with the full covariance (the two estimates are strongly correlated).
# Step 2  contrast    beta_AD ~ a * (b_2h/s_2h - b_fg/s_fg) + c * b_T2D; a > 0 = the AD signal follows
#                     post-load-specific physiology. a equals (th_2h - th_fg)/2 on the same scale only
#                     if b_T2D stands in for the (b_2h/s_2h + b_fg/s_fg) direction.
# Step 3  clusters    Step 2 with the contrast (or b_T2D) interacted with the frozen Suzuki cluster;
#                     headline = Wald heterogeneity test across clusters.
# Estimators: IVW-MRE (weighted least squares, multiplicative random effects floored at 1), which is
# attenuated by noise in the exposure betas, and MRBEE (bias-corrected estimating equations with
# pleiotropy-outlier removal), which corrects for it and for error correlation from sample overlap.
# Step 0 found conditional F(2hGlu) < 10, so Step 2 by MRBEE is PRIMARY and Step 1 is reported as
# underpowered.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# Instruments x outcome on the T2D risk allele, minus the APOE region.
con_outcome_data <- function(inst, g_out, key, apoe) {
  d0 <- as.data.table(inst)[, .(SNP, chr, pos, ea, oa, eaf, beta = beta_t2d, se = se_t2d, p = p_t2d, n = NA_real_)]
  h <- harmonise_pair(d0, g_out, key, "contrast: outcome", exp_name = "t2d", out_name = "ad")
  if (is.null(h$dat)) stop(key, ": no instruments in the outcome GWAS")
  # harmonise_pair() aligns the outcome to the instrument's effect allele (the T2D risk allele).
  y <- as_tibble(h$dat) |> transmute(SNP, beta_y = beta.outcome, se_y = se.outcome)
  d <- inner_join(as_tibble(inst), y, by = "SNP")
  in_apoe <- as.character(d$chr) == as.character(apoe$chr) & d$pos > apoe$start & d$pos < apoe$end
  list(dat = d[!in_apoe, ],
       attrition = bind_rows(h$attrition, attr_row(key, "contrast: outcome", "outside the APOE region",
                                                   sum(!in_apoe), d$SNP[in_apoe])))
}

# Weighted least squares through the origin with a multiplicative random-effects SE (sigma floored
# at 1). Returns coefficients and their covariance.
wls_mre <- function(y, X, se_y) {
  f <- stats::lm(y ~ 0 + X, weights = 1 / se_y^2)
  s <- summary(f)$sigma
  list(b = stats::setNames(stats::coef(f), colnames(X)), V = stats::vcov(f) / s^2 * max(1, s^2), sigma = s,
       n_outlier = NA_integer_)
}

# MRBEE. `R` is the error-correlation matrix ordered (regressors..., outcome).
mrbee_mv <- function(y, X, se_y, seX, R) {
  f <- suppressWarnings(MRBEE::MRBEE.IMRP(by = y, bX = X, byse = se_y, bXse = seX, Rxy = R))
  V <- as.matrix(f$covtheta); dimnames(V) <- list(colnames(X), colnames(X))
  list(b = stats::setNames(as.numeric(f$theta), colnames(X)), V = V, sigma = NA_real_,
       n_outlier = sum(abs(as.numeric(f$delta)) > 1e-8))
}

# Coefficients and linear combinations of a fit as result rows.
fit_rows <- function(fit, combos = list(), n, ...) {
  one <- function(term, w) { b <- sum(w * fit$b); se <- sqrt(drop(t(w) %*% fit$V %*% w))
    tibble(term = term, b = b, se = se, p = 2 * stats::pnorm(-abs(b / se))) }
  k <- length(fit$b)
  bind_rows(c(lapply(seq_len(k), \(i) one(names(fit$b)[i], replace(numeric(k), i, 1))),
              purrr::imap(combos, \(w, nm) one(nm, w)))) |>
    mutate(n_snp = n, ivw_residual_se = fit$sigma, mrbee_n_outlier = fit$n_outlier, ...)
}
# A failed fit is a row with the error message, never a silent drop.
try_fit <- function(expr, step, estimator)
  tryCatch(expr, error = function(e) tibble(term = NA_character_, note = paste("failed:", conditionMessage(e)))) |>
    mutate(step = step, estimator = estimator, .before = 1)

# Error correlations among fasting glucose, 2 h glucose, T2D and the outcome (null-variant z-score
# correlation). MRBEE uses them; non-zero mostly where GWAS share UK Biobank.
con_error_cor <- function(g_fg, g_2h, g_t2d, g_out) {
  g <- list(fg = g_fg, glu2h = g_2h, t2d = g_t2d, y = g_out)
  R <- diag(4); dimnames(R) <- list(names(g), names(g))
  for (i in 1:3) for (j in (i + 1):4) R[i, j] <- R[j, i] <- error_cor(g[[i]], g[[j]])[1, 2]
  R
}

# Steps 1 and 2 on one scale.
con_steps12 <- function(d, R, s_fg, s_2h, scale) {
  bf <- d$beta_fg / s_fg; b2 <- d$beta_2h / s_2h; sf <- d$se_fg / s_fg; s2 <- d$se_2h / s_2h
  n <- nrow(d)
  # Step 1: MVMR
  X1 <- cbind(fg = bf, glu2h = b2); S1 <- cbind(sf, s2); R1 <- R[c("fg", "glu2h", "y"), c("fg", "glu2h", "y")]
  cmb1 <- list(`glu2h - fg` = c(-1, 1))
  # Step 2: contrast + T2D. Error SD of the contrast and its correlation with the other errors follow
  # from the FG/2hGlu error correlation; the per-SNP correlations are averaged (se_2h/se_fg is nearly
  # constant across SNPs, both scaling with 1/sqrt(2pq)).
  xc <- b2 - bf
  sc <- sqrt(s2^2 + sf^2 - 2 * R["fg", "glu2h"] * s2 * sf)
  r_ct <- mean((R["glu2h", "t2d"] * s2 - R["fg", "t2d"] * sf) / sc)
  r_cy <- mean((R["glu2h", "y"] * s2 - R["fg", "y"] * sf) / sc)
  X2 <- cbind(contrast = xc, t2d = d$beta_t2d); S2 <- cbind(sc, d$se_t2d)
  R2 <- matrix(c(1, r_ct, r_cy, r_ct, 1, R["t2d", "y"], r_cy, R["t2d", "y"], 1), 3)
  bind_rows(
    try_fit(fit_rows(wls_mre(d$beta_y, X1, d$se_y), cmb1, n), "1 MVMR", "IVW-MRE"),
    try_fit(fit_rows(mrbee_mv(d$beta_y, X1, d$se_y, S1, R1), cmb1, n), "1 MVMR", "MRBEE"),
    try_fit(fit_rows(wls_mre(d$beta_y, X2, d$se_y), n = n), "2 contrast regression", "IVW-MRE"),
    try_fit(fit_rows(mrbee_mv(d$beta_y, X2, d$se_y, S2, R2), n = n), "2 contrast regression", "MRBEE")) |>
    mutate(scale = scale)
}

# Step 3: per-group slopes of `x` (the contrast, or the T2D effect) with a Wald test of equal slopes.
# Groups with fewer than `min_n` instruments are pooled as "other (small)". IVW-MRE only: MRBEE per
# cluster is too unstable at these instrument counts, so these slopes are attenuated toward zero.
con_groups <- function(d, x, covar, group, label, min_n = 5) {
  g <- as.character(group); g[is.na(g)] <- "unassigned"
  small <- names(which(table(g) < min_n)); g[g %in% small] <- "other (small)"
  lev <- sort(unique(g))
  X <- sapply(lev, \(l) x * (g == l)); colnames(X) <- lev
  Xf <- if (is.null(covar)) X else cbind(X, t2d = covar)
  f <- wls_mre(d$beta_y, Xf, d$se_y)
  k <- length(lev); b <- f$b[1:k]; V <- f$V[1:k, 1:k, drop = FALSE]
  C <- cbind(1, -diag(k - 1)); q <- drop(t(C %*% b) %*% solve(C %*% V %*% t(C)) %*% (C %*% b))
  tibble(model = label, group = lev, n_snp = as.integer(table(g)[lev]), b = b, se = sqrt(diag(V)),
         p = 2 * stats::pnorm(-abs(b / se)), het_Q = q, het_df = k - 1,
         het_p = stats::pchisq(q, k - 1, lower.tail = FALSE))
}

# What the target calls for one instrument set x outcome.
run_contrast_outcome <- function(inst, g_fg, g_2h, g_t2d, g_out, set, outcome, arm, cfg) {
  assert_frozen()
  key <- paste(set, outcome, arm, sep = "__")
  h <- con_outcome_data(inst, g_out, key, cfg$apoe)
  d <- h$dat
  R <- con_error_cor(g_fg, g_2h, g_t2d, g_out)
  sd <- cfg$sd_mmol
  main <- bind_rows(con_steps12(d, R, 1, 1, "per mmol/L"),
                    con_steps12(d, R, sd$fg, sd$glu2h, "per SD")) |>
    mutate(set = set, outcome = outcome, arm = arm, .before = 1)
  a <- readr::read_csv("assignments/instrument_assignments.csv", show_col_types = FALSE)
  d <- left_join(d, a |> select(SNP, cluster), by = "SNP")
  xc <- d$beta_2h / sd$glu2h - d$beta_fg / sd$fg
  beta_cell <- if_else(grepl("^Beta cell", d$cluster), "Beta cell (+PI and -PI)", "all other clusters")
  groups <- if (all(is.na(d$cluster))) NULL else bind_rows(
    con_groups(d, xc, d$beta_t2d, d$cluster, "contrast (per SD) x cluster + T2D"),
    con_groups(d, xc, d$beta_t2d, beta_cell, "contrast (per SD) x beta-cell + T2D"),
    con_groups(d, d$beta_t2d, NULL, d$cluster, "T2D effect x cluster"),
    con_groups(d, d$beta_t2d, NULL, beta_cell, "T2D effect x beta-cell")) |>
    mutate(set = set, outcome = outcome, arm = arm, .before = 1)
  list(rows = main, groups = groups, attrition = h$attrition,
       error_cor = as_tibble(R, rownames = "trait") |> mutate(set = set, outcome = outcome, arm = arm, .before = 1),
       data = d |> mutate(set = set, outcome = outcome, arm = arm, .before = 1))
}

# ---- drug-target layer ----------------------------------------------------------------------------------

# Correlated-instrument IVW (generalised least squares with the LD matrix): log-OR AD per log-OR T2D
# through each target gene and each drug class, plus the FG and 2hGlu effects per log-OR T2D so the
# fasting/post-load profile of each target sits next to its AD estimate.
gls_ratio <- function(bx, by, se_y, Rld) {
  Oi <- solve(outer(se_y, se_y) * Rld)
  v <- 1 / drop(t(bx) %*% Oi %*% bx); b <- v * drop(t(bx) %*% Oi %*% by)
  c(b = b, se = sqrt(v))
}

run_cis_outcome <- function(cis, g_out, set, outcome, arm, cfg, bim, bfile, plink2, min_n = 2) {
  assert_frozen()
  key <- paste("cis", set, outcome, arm, sep = "__")
  d <- con_outcome_data(cis, g_out, key, cfg$apoe)$dat
  one <- function(x, level, name) {
    if (nrow(x) < min_n) return(NULL)
    Rld <- ld_matrix(x$SNP, bim, bfile, plink2)
    x <- x[match(rownames(Rld), x$SNP), ]
    s <- if_else(x$ea == bim$a1[match(x$SNP, bim$SNP)], 1, -1)      # ld_matrix() is signed on bim allele a1
    Rld <- Rld * outer(s, s)
    est <- function(by, se, trait) { r <- gls_ratio(x$beta_t2d, by, se, Rld)
      tibble(trait = trait, b = r[["b"]], se = r[["se"]], p = 2 * stats::pnorm(-abs(b / se))) }
    bind_rows(est(x$beta_y, x$se_y, "AD (log-OR)"), est(x$beta_fg, x$se_fg, "fasting glucose (mmol/L)"),
              est(x$beta_2h, x$se_2h, "2h glucose (mmol/L)")) |>
      mutate(level = level, target = name, n_snp = nrow(x), max_abs_r = max(abs(Rld[upper.tri(Rld)])),
             within_2Mb_of_APOE = any(as.character(x$chr) == "19" & abs(x$pos - 45.41e6) < 2e6), .before = 1)
  }
  bind_rows(unname(purrr::compact(c(
    purrr::imap(split(d, d$drug_target_gene), \(x, nm) one(x, "gene", nm)),
    purrr::imap(split(d, d$drug_class), \(x, nm) one(x, "class", nm)))))) |>
    mutate(set = set, outcome = outcome, arm = arm, scale = "per log-odds of T2D", .before = 1)
}
