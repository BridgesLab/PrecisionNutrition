#!/usr/bin/env Rscript
# run_cause_greatlakes.R --------------------------------------------------------
# Headless CAUSE fit for one exposure arm. Does exactly what the merge chunk in
# robust_mr_cause.qmd does, but without knitr, so it can be sbatch'd.
#
#   Rscript scripts/run_cause_greatlakes.R <arm> [project_root]
#
#   <arm>  one of: overlapping | overlap-free
#
# Writes: <cache>/cause_fit_<arm>.rds  — the qmd picks this up automatically.

suppressPackageStartupMessages({
  library(tidyverse)
  library(cause)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_cause_greatlakes.R <arm> [project_root]")
arm  <- args[1]
root <- if (length(args) >= 2) normalizePath(args[2]) else getwd()

stopifnot(arm %in% c("overlapping", "overlap-free"))

source(file.path(root, "R", "robust_mr_helpers.R"))
set_project_root(root)
cfg <- load_robust_cfg(file.path(root, "config_robust_mr.yml"))
set.seed(cfg$seed)

exp_name <- if (arm == "overlapping") "glgc2021_ldl" else "willer2013_ldl"
cache    <- file.path(root, cfg$paths$cache)
fit_path <- file.path(cache, paste0("cause_fit_", arm, ".rds"))

if (file.exists(fit_path)) {
  message("Fit already exists at ", fit_path, " — nothing to do.")
  quit(save = "no", status = 0)
}

t0 <- Sys.time()
msg <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

msg("Reading harmonised inputs")
e <- readRDS(file.path(cache, paste0(exp_name, "_cause.rds")))
o <- readRDS(file.path(cache, "ebmd_cause.rds"))
msg(sprintf("  exposure: %s SNPs | outcome: %s SNPs",
            format(nrow(e), big.mark = ","), format(nrow(o), big.mark = ",")))

msg("gwas_merge() — this is the memory peak")
X <- cause::gwas_merge(
  e, o,
  snp_name_cols = c("SNP", "SNP"),
  beta_hat_cols = c("beta_hat", "beta_hat"),
  se_cols       = c("se", "se"),
  A1_cols       = c("A1", "A1"),
  A2_cols       = c("A2", "A2")
)
rm(e, o); gc()
msg(sprintf("  merged: %s variants", format(nrow(X), big.mark = ",")))

msg("est_cause_params() on ", format(cfg$cause$n_param_variants, big.mark = ","), " variants")
varlist <- with(X, sample(snp, size = min(cfg$cause$n_param_variants, nrow(X)),
                          replace = FALSE))
params  <- cause::est_cause_params(X, varlist)
msg(sprintf("  rho = %.4f", params$rho))

msg("LD pruning at p < ", cfg$cause$prune_p, ", r2 < ", cfg$cause$prune_r2)
X <- X %>% dplyr::mutate(pval1 = 2 * pnorm(abs(beta_hat_1 / seb1), lower.tail = FALSE))
top <- clump_local(
  X, snp_col = "snp", p_col = "pval1",
  r2 = cfg$cause$prune_r2, kb = cfg$cause$prune_kb, p_thresh = cfg$cause$prune_p,
  bfile     = file.path(root, cfg$paths$plink_bfile),
  plink_bin = cfg$paths$plink_bin
)
msg(sprintf("  %s variants retained", format(nrow(top), big.mark = ",")))

msg("cause() — fitting null / sharing / causal")
fit <- cause::cause(X = X, variants = top$snp, param_ests = params)

dir.create(cache, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(fit = fit, params = params,
             n_variants = nrow(top), n_merged = nrow(X)),
        fit_path)

msg("Done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
msg("Wrote ", fit_path)
print(fit$elpd)
