# latent.R — latent "residual heritable confounding" terms from CAUSE and MR-APSS.
#
# Spec: one latent-confounding term per association and method, with an interval, next to the causal
# estimate. It covers ALL heritable shared factors (not a named confounder; those go through MVMR
# triads); purely environmental confounding does not bias MR and is not in it.
#
# Exact object paths used (read from the package sources, 2026-10-07; cause 1.2.0, MRAPSS 0.0.0.9000):
#   CAUSE  fit$causal$joint_post   grid posterior: gamma, eta, q, <p>start, <p>stop, log_post
#          fit$sharing$joint_post  grid posterior: eta, q (gamma fixed at 0)
#          cause:::samp_from_grid(joint_post, params, n)   CAUSE's own joint sampler (cell by
#                                  posterior weight, then uniform within the cell)
#          fit$elpd                row model1 == "sharing", model2 == "causal": delta_elpd, se_delta_elpd
#   MRAPSS est_paras()$Omega       = ldsc_GC()$cov / M   (background covariance, per-SNP, on the
#                                  z/sqrt(N) scale that b.exp/b.out and beta use)
#          est_paras()$C           = ldsc_GC()$I         (LDSC intercepts: sample structure)
#          inside ldsc_GC():       h2_1$delete.values[, 1], rho_g$delete.values[, 1], n.blocks (200)
#                                  -- computed by the package but not returned; est_paras_jk() below
#                                  retains them (no other change).

suppressPackageStartupMessages({ library(dplyr) })

# ---- MR-APSS: est_paras() with the LDSC jackknife delete-values retained -------------------------
# ldsc_GC() is run unchanged except its final return(), which also returns the per-block delete-one
# estimates; est_paras() is run in an environment where that version is the ldsc_GC() it finds.
# check_est_paras_jk() confirms C and Omega are identical to the unmodified package.
ldsc_GC_jk <- local({
  f <- get("ldsc_GC", envir = asNamespace("MRAPSS"))
  b <- body(f); n <- length(b)
  stopifnot(identical(b[[n]][[1]], as.name("return")))
  b[[n]] <- quote(return(list(cov = cov, cov.se = cov.se, I = I, I.se = I.se, rg = rg, rg.se = rg.se,
                              jk = list(n.blocks = n.blocks,
                                        h2_1 = h2_1$delete.values[, 1], gcov = rho_g$delete.values[, 1]))))
  body(f) <- b
  environment(f) <- asNamespace("MRAPSS")
  f
})

est_paras_jk <- local({
  f <- get("est_paras", envir = asNamespace("MRAPSS"))
  env <- new.env(parent = asNamespace("MRAPSS"))
  assign("ldsc_GC", ldsc_GC_jk, envir = env)
  environment(f) <- env
  f
})

# Background slope Omega_xy / Omega_xx with a block-jackknife interval. Depending on the LDSC path,
# the delete-one values are either already on the cov scale (the two-step path used by est_paras:
# checked, delete-one means 0.0575 / 0.0112 vs cov 0.0574 / 0.0113) or are raw regression slopes
# (cov = slope / N.bar * M). Each series is therefore rescaled by full / mean(delete-one), which is
# ~1 in the first case and the constant M / N.bar in the second; the ratio's pseudo-values follow.
# (The package's own rg.se uses the delete-one values unscaled.)
jk_rescale <- function(del, full) del * full / mean(del)

apss_background_slope <- function(paras, n_exp, n_out) {
  O <- paras$Omega; C <- paras$C; jk <- paras$ldsc_res$jk
  slope <- O[1, 2] / O[1, 1]
  cov <- paras$ldsc_res$cov
  ratio_del <- jk_rescale(jk$gcov, cov[1, 2]) / jk_rescale(jk$h2_1, cov[1, 1])
  B <- jk$n.blocks
  pseudo <- B * slope - (B - 1) * ratio_del
  se <- sqrt(stats::var(pseudo) / B)
  tibble(apss_bg_slope = slope, apss_bg_slope_se = se,
         apss_bg_slope_lo = slope - 1.96 * se, apss_bg_slope_hi = slope + 1.96 * se,
         apss_omega11 = O[1, 1], apss_omega12 = O[1, 2], apss_omega22 = O[2, 2],
         # Sample-structure slope: C_xy / C_xx on the z scale, and rescaled to the beta scale.
         # b = z / sqrt(N), so noise cov is C12 / sqrt(N1 N2) and exposure noise var is C11 / N1:
         # slope on the beta scale = (C12 / C11) * sqrt(N1 / N2).
         apss_C11 = C[1, 1], apss_C12 = C[1, 2], apss_C22 = C[2, 2],
         apss_C_slope_z = C[1, 2] / C[1, 1],
         apss_C_slope_beta = (C[1, 2] / C[1, 1]) * sqrt(n_exp / n_out))
}

# One-off check (run on any pair): the retained-jackknife version must reproduce the package.
check_est_paras_jk <- function(e, o, ldsc) {
  a <- MRAPSS::est_paras(dat1 = e, dat2 = o, ld = ldsc$ld, M = ldsc$M)
  b <- est_paras_jk(dat1 = e, dat2 = o, ld = ldsc$ld, M = ldsc$M)
  # Rescaled delete-one cov12 values must reproduce the package's own jackknife SE for cov12.
  jk <- b$ldsc_res$jk; B <- jk$n.blocks
  g <- jk_rescale(jk$gcov, b$ldsc_res$cov[1, 2])
  pseudo <- B * b$ldsc_res$cov[1, 2] - (B - 1) * g
  tibble(C_identical = isTRUE(all.equal(a$C, b$C)), Omega_identical = isTRUE(all.equal(a$Omega, b$Omega)),
         gcov_se_package = b$ldsc_res$cov.se[1, 2], gcov_se_from_deletes = sqrt(stats::var(pseudo) / B))
}

# ---- CAUSE: joint posterior of (gamma, eta, q) ----------------------------------------------------
summ <- function(x, pre) {
  q <- stats::quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE)
  setNames(as.list(unname(q)), paste0(pre, c("_est", "_lo", "_hi")))
}

# Latent confounding slope q*eta and the fraction of the naive slope it accounts for,
# q*eta / (gamma + q*eta), computed per JOINT draw under the causal model.
cause_latent <- function(post_causal, post_sharing, elpd, nsamps = 20000, seed = 1) {
  set.seed(seed)
  d <- cause:::samp_from_grid(post_causal, c("gamma", "eta", "q"), nsamps)
  conf <- d$q * d$eta
  denom <- d$gamma + conf
  frac <- conf / denom
  s <- cause:::samp_from_grid(post_sharing, c("eta", "q"), nsamps)
  row <- if (!is.null(elpd)) elpd[elpd$model1 == "sharing" & elpd$model2 == "causal", ] else NULL

  out <- as_tibble(c(summ(d$gamma, "cause_gamma"), summ(d$eta, "cause_eta"), summ(d$q, "cause_q"),
                     summ(conf, "cause_qeta"), summ(denom, "cause_naive"), summ(frac, "cause_frac"),
                     summ(s$q, "cause_q_sharing")))
  denom_spans0 <- out$cause_naive_lo < 0 & out$cause_naive_hi > 0
  q_high <- out$cause_q_est > 0.5
  out |> mutate(
    cause_delta_elpd = if (!is.null(row) && nrow(row)) row$delta_elpd[1] else NA_real_,
    cause_delta_elpd_se = if (!is.null(row) && nrow(row)) row$se_delta_elpd[1] else NA_real_,
    cause_flags = paste(c(if (q_high) "q > 0.5: causal and sharing poorly separated, fraction unreliable",
                          if (denom_spans0) "gamma + q*eta interval spans 0: fraction not reported"),
                        collapse = "; "),
    # Spec: suppress the fraction when its denominator's interval spans zero.
    across(c(cause_frac_est, cause_frac_lo, cause_frac_hi), \(x) if_else(denom_spans0, NA_real_, x)))
}

# ---- one row per association x method, in the spec's column layout --------------------------------
latent_table <- function(cause_rows, apss_rows) {
  # Either input may be empty (e.g. locally, where CAUSE cannot run); keep whatever is present.
  has <- function(d, col) !is.null(d) && nrow(d) > 0 && col %in% names(d)
  cz <- if (!has(cause_rows, "cause_qeta_est")) NULL else cause_rows |>
    transmute(pair_id = assoc_id, method = "CAUSE",
              causal_est = cause_gamma_est, causal_lo = cause_gamma_lo, causal_hi = cause_gamma_hi,
              latent_conf_est = cause_qeta_est, latent_conf_lo = cause_qeta_lo, latent_conf_hi = cause_qeta_hi,
              latent_conf_label = "q*eta: latent shared-factor slope (CAUSE causal model, joint posterior)",
              frac_shared_est = cause_frac_est, frac_shared_lo = cause_frac_lo, frac_shared_hi = cause_frac_hi,
              cause_q = cause_q_est, cause_q_sharing = cause_q_sharing_est,
              cause_delta_elpd, cause_delta_elpd_se,
              apss_sample_structure_slope = NA_real_, flags = cause_flags)
  ap <- if (!has(apss_rows, "apss_bg_slope")) NULL else apss_rows |>
    transmute(pair_id = assoc_id, method = "MR-APSS",
              causal_est = b, causal_lo = ci_lo, causal_hi = ci_hi,
              latent_conf_est = apss_bg_slope, latent_conf_lo = apss_bg_slope_lo, latent_conf_hi = apss_bg_slope_hi,
              latent_conf_label = paste("background (non-foreground) slope Omega_xy/Omega_xx: includes residual",
                                        "heritable confounding; may include some causal signal from non-IV SNPs"),
              frac_shared_est = NA_real_, frac_shared_lo = NA_real_, frac_shared_hi = NA_real_,
              cause_q = NA_real_, cause_q_sharing = NA_real_, cause_delta_elpd = NA_real_,
              cause_delta_elpd_se = NA_real_,
              apss_sample_structure_slope = apss_C_slope_beta,
              flags = if_else(apss_bg_slope_lo < 0 & apss_bg_slope_hi > 0, "background slope interval spans 0", ""))
  out <- bind_rows(cz, ap)
  if (nrow(out)) arrange(out, pair_id, method) else out
}
