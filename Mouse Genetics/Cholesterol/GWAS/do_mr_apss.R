## MR-APSS: bidirectional MR for DO mouse cholesterol ↔ calcium
##
## Run on Great Lakes after both GEMMA GWAS are complete.
## Writes do_mr_apss_results.csv; download to local GWAS/ folder then re-run
## do_mr_mrbee.R to produce the combined table.
##
## Install dependencies (once):
##   module load R/4.3.1
##   Rscript -e "remotes::install_github('stephenslab/mixsqp')"
##   Rscript -e "remotes::install_github('stephens999/ashr')"
##   Rscript -e "remotes::install_github('YangLabHKUST/MR-APSS')"
##
## C/Omega approximation rationale:
##   MR-APSS normally estimates C (LDSC intercept matrix) and Omega (genetic
##   covariance) from bivariate LD score regression, which requires human LD
##   scores.  DO mice lack a published LD score reference, so both are estimated
##   directly from the GWAS summary stats:
##     C[i,i]  = lambda_GC_i  (genomic inflation; captures relatedness/overlap)
##     C[1,2]  = mean cross-trait Z product for null SNPs (the LDSC cross-trait
##               intercept concept, same quantity MRBEE::errorCov uses)
##     Omega   = mean chi2 excess above C, normalised by N
##   L2 is set to 1 for all instruments (DO LD diversity ≈ uniform LD scores).

suppressPackageStartupMessages({
  library(tidyverse)
  library(MRAPSS)    # YangLabHKUST/MR-APSS
})

# ── Paths ────────────────────────────────────────────────────────────────────
GEMMA_DIR <- "/home/davebrid/bridgeslab/GEMMA/DO_HFD"
OUT_DIR   <- file.path(GEMMA_DIR, "output")
N_CA      <- 763L   # non-NA calcium2 sample size (from extract_calcium_phenotype.R)

# ── Thresholds ───────────────────────────────────────────────────────────────
P_THRESH <- 5e-6
WIND_MB  <- 1
MIN_F    <- 10

# ── Load and harmonise GWAS summary stats ────────────────────────────────────
read_gemma <- function(path, trait) {
  read_tsv(path, col_types = cols(.default = "d", rs = "c",
                                  allele1 = "c", allele0 = "c"),
           show_col_types = FALSE) |>
    rename(snp = rs, pos = ps, a1 = allele1, a0 = allele0, pval = p_wald) |>
    mutate(trait = trait, chr_num = as.integer(str_extract(snp, "^[^_]+"))) |>
    filter(!is.na(beta), !is.na(se), se > 0)
}

chol <- read_gemma(file.path(OUT_DIR, "cholesterol_all.assoc.txt"), "cholesterol")
ca   <- read_gemma(file.path(OUT_DIR, "calcium2_all.assoc.txt"),    "calcium2")

harm <- inner_join(
  chol |> select(snp, chr_num, pos, a1, a0,
                 beta_chol = beta, se_chol = se, pval_chol = pval),
  ca   |> select(snp, beta_ca = beta, se_ca = se, pval_ca = pval),
  by = "snp"
) |>
  mutate(F_chol = (beta_chol / se_chol)^2,
         F_ca   = (beta_ca   / se_ca)^2,
         z_chol = beta_chol / se_chol,
         z_ca   = beta_ca   / se_ca)

cat("Harmonised SNPs:", nrow(harm), "\n")

# ── Estimate C matrix (no DO mouse LD scores available) ──────────────────────
lambda_chol <- median(harm$z_chol^2) / 0.4549
lambda_ca   <- median(harm$z_ca^2)   / 0.4549

null <- harm |> filter(pval_chol > 0.05, pval_ca > 0.05)
c12  <- mean(null$z_chol * null$z_ca)   # cross-trait LDSC-intercept equivalent

C <- matrix(c(lambda_chol, c12, c12, lambda_ca), nrow = 2)
cat("C matrix (genomic inflation + cross-trait intercept):\n")
print(round(C, 4))
if (!all(eigen(C)$values > 0)) {
  warning("C is not positive definite; setting off-diagonal to 0")
  C[1, 2] <- C[2, 1] <- 0
}

# ── Estimate Omega matrix ─────────────────────────────────────────────────────
# Mean chi2 excess above C, normalised by N (with L2=1 this is the h2 per SNP)
omega_11 <- max(1e-6, (mean(harm$z_chol^2) - lambda_chol) / N_CA)
omega_22 <- max(1e-6, (mean(harm$z_ca^2)   - lambda_ca)   / N_CA)
omega_12 <- (mean(harm$z_chol * harm$z_ca) - c12) / N_CA

Omega <- matrix(c(omega_11, omega_12, omega_12, omega_22), nrow = 2)
if (!all(eigen(Omega)$values > 0)) {
  message("Omega not PD; using diagonal approximation")
  Omega <- diag(c(omega_11, omega_22))
}
cat("Omega matrix:\n"); print(signif(Omega, 3))

# ── Instrument selection: top SNP per Mb window ──────────────────────────────
select_instruments <- function(df, pval_col, f_col,
                               p_thresh = P_THRESH, wind_mb = WIND_MB,
                               min_f = MIN_F) {
  sig <- df |> filter(.data[[pval_col]] < p_thresh) |> arrange(.data[[pval_col]])
  if (nrow(sig) == 0) {
    warning("No SNPs below p < ", p_thresh, " for ", pval_col)
    return(sig[0, ])
  }
  kept <- list()
  while (nrow(sig) > 0) {
    top <- sig[1, ]
    kept[[length(kept) + 1]] <- top
    sig <- sig |> filter(!(chr_num == top$chr_num & abs(pos - top$pos) < wind_mb * 1e6))
  }
  bind_rows(kept) |> filter(.data[[f_col]] >= min_f)
}

iv_chol <- select_instruments(harm, "pval_chol", "F_chol")
iv_ca   <- select_instruments(harm, "pval_ca",   "F_ca")
cat("Cholesterol instruments:", nrow(iv_chol), "\n")
cat("Calcium instruments:",     nrow(iv_ca),   "\n")

# ── Build MRdat for MRAPSS() ─────────────────────────────────────────────────
# MR-APSS expects: SNP, A1, A2, b.exp, b.out, se.exp, se.out,
#                  pval.exp, pval.out, L2, Threshold, chi2.exp, chi2.out
# L2 = 1 (uniform LD score approximation for DO mice)
make_mrdat <- function(ivs, exp_b, exp_se, out_b, out_se, pval_exp, pval_out,
                       thresh) {
  ivs |>
    transmute(SNP      = snp,
              A1       = a1,
              A2       = a0,
              b.exp    = .data[[exp_b]],
              se.exp   = .data[[exp_se]],
              b.out    = .data[[out_b]],
              se.out   = .data[[out_se]],
              pval.exp = .data[[pval_exp]],
              pval.out = .data[[pval_out]],
              chi2.exp = (.data[[exp_b]] / .data[[exp_se]])^2,
              chi2.out = (.data[[out_b]] / .data[[out_se]])^2,
              L2        = 1.0,
              Threshold = as.numeric(thresh)) |>
    filter(complete.cases(pick(b.exp, b.out, se.exp, se.out, L2))) |>
    as.data.frame()
}

# ── Run MR-APSS (one direction) ──────────────────────────────────────────────
run_apss <- function(mrdat, C, Omega, label) {
  if (nrow(mrdat) < 4) {
    warning(label, ": fewer than 4 instruments, skipping MR-APSS")
    return(NULL)
  }

  res <- tryCatch(
    MRAPSS::MRAPSS(mrdat, exposure = "exposure", outcome = "outcome",
                   C = C, Omega = Omega),
    error = function(e) { warning(label, " MR-APSS error: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)

  pluck_num <- function(x, ...) {
    for (nm in c(...)) {
      v <- x[[nm]]
      if (!is.null(v)) { out <- suppressWarnings(as.numeric(v[1])); if (!is.na(out)) return(out) }
    }
    NA_real_
  }

  b  <- pluck_num(res, "beta")
  se <- pluck_num(res, "beta.se")
  tibble(arm       = label,
         method    = "MR-APSS",
         n_ivs     = nrow(mrdat),
         mean_F    = round(mean((mrdat$b.exp / mrdat$se.exp)^2), 1),
         beta      = b,
         se        = se,
         ci_lo     = b - 1.96 * se,
         ci_hi     = b + 1.96 * se,
         pval      = suppressWarnings(as.numeric(res$pvalue)),
         n_pleiotropic = if (!is.null(res$post$Pi))
                           sum(as.numeric(res$post$Pi) < 0.5) else NA_integer_,
         pi0       = pluck_num(res, "pi0"),
         tau2      = pluck_num(res, "tau.sq"))
}

# ── Bidirectional ─────────────────────────────────────────────────────────────
cat("\n--- Cholesterol → Calcium ---\n")
mrdat_chol <- make_mrdat(iv_chol,
                         "beta_chol", "se_chol", "beta_ca", "se_ca",
                         "pval_chol", "pval_ca", P_THRESH)
r1 <- run_apss(mrdat_chol, C, Omega, "cholesterol → calcium")

cat("\n--- Calcium → Cholesterol ---\n")
mrdat_ca <- make_mrdat(iv_ca,
                       "beta_ca", "se_ca", "beta_chol", "se_chol",
                       "pval_ca", "pval_chol", P_THRESH)
r2 <- run_apss(mrdat_ca, C, Omega, "calcium → cholesterol")

results <- bind_rows(r1, r2)
cat("\n=== MR-APSS results ===\n")
print(results, width = Inf)

write_csv(results, file.path(GEMMA_DIR, "do_mr_apss_results.csv"))
cat("\nWritten → do_mr_apss_results.csv\n")
cat("Download to local GWAS/ folder then re-run do_mr_mrbee.R to combine.\n")
