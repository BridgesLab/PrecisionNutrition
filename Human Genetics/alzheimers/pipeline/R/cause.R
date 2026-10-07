# cause.R — CAUSE (Morrison et al. 2020) for associations whose GWAS share samples.
#
# Follows calcium-cholesterol/scripts/run_cause_greatlakes.R, which ran successfully on Great Lakes:
# HapMap3-restricted input on a per-SD scale, nuisance parameters on up to 1M random variants, LD
# pruning at p < 1e-3 / r2 < 0.01, then the null/sharing/causal fits. The posterior summaries reuse
# calcium-cholesterol/R/robust_mr_helpers.R (tidy_cause, cause_posterior_quantiles, cause_p_below,
# cause_model_weights_boot), sourced by _targets.R.
#
# Decision rule (collaborator's recommendation): the 95% credible interval for gamma and
# P(gamma < 0) under the causal model are PRIMARY; the sharing-vs-causal ELPD test is reported
# alongside. CAUSE needs loo < 2.10 (cause 1.2.0 indexes loo_compare() positionally), so these
# targets run on Great Lakes, where loo 2.4.1 is pinned.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# Standard format -> CAUSE input, per-SD scale (z / sqrt(N)), HapMap3 only.
to_cause_input <- function(g, hm3, n_fallback) {
  g[SNP %chin% hm3][, nn := fifelse(is.na(n), n_fallback, n)][
    , .(SNP, A1 = ea, A2 = oa, beta_hat = (beta / se) / sqrt(nn), se = 1 / sqrt(nn), pval = p)]
}

run_cause <- function(g_exp, g_out, key, hm3, n_exp, n_out, bfile, plink2, cfg, seed = 20261004) {
  if (packageVersion("loo") >= "2.10")
    stop("CAUSE needs loo < 2.10 (installed: ", packageVersion("loo"), "). Run on Great Lakes, ",
         "or: remotes::install_version('loo', '2.4.1')")
  set.seed(seed)
  e <- to_cause_input(g_exp, hm3, n_exp); o <- to_cause_input(g_out, hm3, n_out)
  X <- cause::gwas_merge(e, o, snp_name_cols = c("SNP", "SNP"),
                         beta_hat_cols = c("beta_hat", "beta_hat"), se_cols = c("se", "se"),
                         A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))
  att <- attr_row(key, "cause", "HapMap3, merged with outcome (gwas_merge drops ambiguous/mismatched)",
                  nrow(X))
  vars <- sample(X$snp, min(cfg$n_param_variants, nrow(X)))
  params <- cause::est_cause_params(X, vars)

  X$pval1 <- 2 * stats::pnorm(-abs(X$beta_hat_1 / X$seb1))
  cand <- data.table(SNP = X$snp, p = X$pval1)[p < cfg$prune_p]
  top <- plink_clump(cand, bfile, cfg$prune_p, cfg$prune_r2, cfg$prune_kb, plink2)
  att <- bind_rows(att,
    attr_row(key, "cause", paste0("p < ", cfg$prune_p), nrow(cand)),
    attr_row(key, "cause", sprintf("LD prune r2 < %s, %s kb", cfg$prune_r2, cfg$prune_kb),
             length(top), setdiff(cand$SNP, top)))

  fit <- cause::cause(X = X, variants = top, param_ests = params)
  qtab <- cause_posterior_quantiles(fit, "causal")
  td <- tidy_cause(fit, exposure = key, outcome = key, arm = key, ci_size = 0.95)
  res <- tibble(method = "CAUSE", scale = "SD(outcome) per SD(exposure)",
                b = td$b, ci_lo = td$lci, ci_hi = td$uci,
                p_gamma_below_0 = cause_p_below(qtab, "gamma", 0),
                ci_excludes_0 = td$lci > 0 | td$uci < 0,
                n_snp = length(top), cause_rho = params$rho,
                cause_eta_median = td$eta_med, cause_q_median = td$q_med,
                elpd_sharing_vs_causal = td$delta_elpd, elpd_se = td$se_delta_elpd,
                elpd_p = td$pval)
  # Latent confounding (pipeline/R/latent.R): q*eta and q*eta/(gamma + q*eta) per joint posterior
  # draw. The grid posteriors are kept so this never needs another CAUSE run.
  post <- list(causal = fit$causal$joint_post, sharing = fit$sharing$joint_post)
  res <- bind_cols(res, cause_latent(post$causal, post$sharing, as.data.frame(fit$elpd), seed = seed))
  inst <- as.data.table(X)[snp %chin% top, .(SNP = snp, effect_allele = A1, other_allele = A2,
                                             beta_exposure = beta_hat_1, se_exposure = seb1,
                                             p_exposure = pval1, beta_outcome = beta_hat_2,
                                             se_outcome = seb2)]
  list(result = res, instruments = inst, attrition = att,
       weights = cause_model_weights_boot(as.data.frame(fit$elpd)), elpd = as.data.frame(fit$elpd),
       posterior = post)
}
