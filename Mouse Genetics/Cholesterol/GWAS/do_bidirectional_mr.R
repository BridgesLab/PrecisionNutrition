## Bidirectional MR: cholesterol ↔ calcium in DO mice (Svenson HFD v12)
##
## Same individuals provide both traits → 100% sample overlap.
## Uses MRBEE (bias-corrected estimator with overlap correction) as primary.
## IVW is included for reference but will be biased; treat as directional check.
##
## Run on Great Lakes after both GWAS are complete:
##   module load R/4.3.1
##   Rscript do_bidirectional_mr.R
##
## Expects in GEMMA_DIR/output/:
##   cholesterol_all.assoc.txt   (cholesterol GWAS)
##   calcium2_all.assoc.txt      (19-week calcium GWAS, primary)

library(tidyverse)
library(MRBEE)       # install: devtools::install_github("siqixu/MRBEE")

# ── Paths ────────────────────────────────────────────────────────────────────
GEMMA_DIR <- "/home/davebrid/bridgeslab/GEMMA/DO_HFD"
OUT_DIR   <- file.path(GEMMA_DIR, "output")

# ── Thresholds ───────────────────────────────────────────────────────────────
P_THRESH  <- 5e-6     # instrument selection (genome-wide for mice ~5e-6)
WIND_MB   <- 1        # one SNP per Mb window for approximate LD pruning
MIN_F     <- 10       # minimum F-statistic for instrument strength check

# ── Load GWAS summary stats ───────────────────────────────────────────────────
read_gemma <- function(path, trait) {
  read_tsv(path, col_types = cols(.default = "d", rs = "c",
                                  allele1 = "c", allele0 = "c")) |>
    rename(snp = rs, pos = ps, a1 = allele1, a0 = allele0,
           beta = beta, se = se, pval = p_wald) |>
    mutate(trait = trait,
           # chr is -9 in GEMMA output when annotation is missing;
           # parse from SNP id (format: chr_pos_a0_a1)
           chr_num = as.integer(str_extract(snp, "^[^_]+"))) |>
    filter(!is.na(beta), !is.na(se), se > 0)
}

chol <- read_gemma(file.path(OUT_DIR, "cholesterol_all.assoc.txt"), "cholesterol")
ca   <- read_gemma(file.path(OUT_DIR, "calcium2_all.assoc.txt"),    "calcium2")

cat("Cholesterol GWAS SNPs:", nrow(chol), "\n")
cat("Calcium GWAS SNPs:", nrow(ca), "\n")

# ── Harmonise: inner join on SNP id ──────────────────────────────────────────
# Both traits are rank-normalised from the same sample, so alleles are identical
harm <- inner_join(
  chol |> select(snp, chr_num, pos, a1, a0, af,
                 beta_chol = beta, se_chol = se, pval_chol = pval),
  ca   |> select(snp, beta_ca = beta, se_ca = se, pval_ca = pval),
  by = "snp"
) |>
  # F-statistic approximation: (beta/se)^2 for single-SNP strength
  mutate(F_chol = (beta_chol / se_chol)^2,
         F_ca   = (beta_ca   / se_ca)^2)

cat("Harmonised SNPs:", nrow(harm), "\n")

# ── Instrument selection: top SNP per Mb window ──────────────────────────────
select_instruments <- function(df, pval_col, beta_col, se_col,
                               p_thresh = P_THRESH, wind_mb = WIND_MB,
                               min_f    = MIN_F) {
  sig <- df |>
    filter(.data[[pval_col]] < p_thresh) |>
    arrange(.data[[pval_col]])

  if (nrow(sig) == 0) {
    warning("No SNPs below p < ", p_thresh, " for ", pval_col)
    return(sig)
  }

  # Greedy window clump: keep top SNP, remove within wind_mb, repeat
  kept <- list()
  while (nrow(sig) > 0) {
    top <- sig[1, ]
    kept[[length(kept) + 1]] <- top
    sig <- sig |>
      filter(!(chr_num == top$chr_num &
                 abs(pos - top$pos) < wind_mb * 1e6))
  }
  ivs <- bind_rows(kept) |>
    filter(.data[[ sub("pval_", "F_", pval_col) ]] >= min_f)
  cat("Instruments (", pval_col, "): ", nrow(ivs), " SNPs\n", sep = "")
  ivs
}

iv_chol <- select_instruments(harm, "pval_chol", "beta_chol", "se_chol")
iv_ca   <- select_instruments(harm, "pval_ca",   "beta_ca",   "se_ca")

# ── Phenotypic correlation (same-sample overlap correction for MRBEE) ─────────
# Load raw phenotype files to compute observed r(cholesterol, calcium)
chol_pheno <- as.numeric(readLines(file.path(GEMMA_DIR, "Cholesterol_all.txt")))
ca_pheno   <- as.numeric(readLines(file.path(GEMMA_DIR, "Calcium2_all.txt")))
# Both are rank-normalised; drop NAs before correlating
ok <- !is.na(chol_pheno) & !is.na(ca_pheno)
rho <- cor(chol_pheno[ok], ca_pheno[ok], method = "pearson")
cat("Phenotypic correlation (chol, ca):", round(rho, 4), "\n")

# MRBEE Sigma_e (2×2 residual covariance; diagonal = 1 for normalised traits)
Sigma_e <- matrix(c(1, rho, rho, 1), nrow = 2)

# ── MRBEE helper ─────────────────────────────────────────────────────────────
run_mrbee <- function(ivs, exp_beta, exp_se, out_beta, out_se, label) {
  if (nrow(ivs) < 3) {
    warning(label, ": fewer than 3 instruments, skipping MRBEE")
    return(NULL)
  }

  bx <- ivs[[exp_beta]]
  bxse <- ivs[[exp_se]]
  by <- ivs[[out_beta]]
  byse <- ivs[[out_se]]

  # Instrument strength summary
  F_mean <- mean((bx / bxse)^2)
  cat(label, "— n instruments:", nrow(ivs), "| mean F:", round(F_mean, 1), "\n")

  # MRBEE.IPOD: iterative pleiotropy-robust estimator with overlap correction
  res <- tryCatch(
    MRBEE.IPOD(
      by  = by,
      bx  = matrix(bx, ncol = 1),
      byse = byse,
      bxse = matrix(bxse, ncol = 1),
      Sigma_e = Sigma_e,
      tauvec  = seq(0, 0.5, by = 0.05)   # grid for pleiotropy variance
    ),
    error = function(e) { warning(label, " MRBEE.IPOD failed: ", e$message); NULL }
  )

  if (is.null(res)) return(NULL)

  tibble(
    arm        = label,
    method     = "MRBEE.IPOD",
    n_ivs      = nrow(ivs),
    mean_F     = round(F_mean, 1),
    beta       = res$theta,
    se         = res$se,
    pval       = 2 * pnorm(abs(res$theta / res$se), lower.tail = FALSE),
    n_outliers = res$nout
  )
}

# IVW (naive, ignores overlap — reference only)
run_ivw <- function(ivs, exp_beta, exp_se, out_beta, out_se, label) {
  if (nrow(ivs) < 2) return(NULL)
  bx   <- ivs[[exp_beta]];  bxse <- ivs[[exp_se]]
  by   <- ivs[[out_beta]];  byse <- ivs[[out_se]]
  w    <- 1 / byse^2
  beta <- sum(w * by * bx / bx^2) / sum(w / bx^2)   # IVW fixed-effects
  # corrected SE via delta method
  se   <- sqrt(1 / sum(w / bx^2))
  tibble(arm = label, method = "IVW (naive)",
         n_ivs = nrow(ivs), mean_F = round(mean((bx/bxse)^2), 1),
         beta = beta, se = se,
         pval = 2 * pnorm(abs(beta/se), lower.tail = FALSE),
         n_outliers = NA_integer_)
}

# ── Run bidirectional MR ──────────────────────────────────────────────────────
results <- bind_rows(
  run_mrbee(iv_chol, "beta_chol", "se_chol", "beta_ca", "se_ca",
            "cholesterol → calcium"),
  run_ivw(  iv_chol, "beta_chol", "se_chol", "beta_ca", "se_ca",
            "cholesterol → calcium"),
  run_mrbee(iv_ca,   "beta_ca", "se_ca", "beta_chol", "se_chol",
            "calcium → cholesterol"),
  run_ivw(  iv_ca,   "beta_ca", "se_ca", "beta_chol", "se_chol",
            "calcium → cholesterol")
)

results <- results |>
  mutate(ci_lo = beta - 1.96 * se,
         ci_hi = beta + 1.96 * se) |>
  relocate(arm, method, n_ivs, mean_F, beta, se, ci_lo, ci_hi, pval, n_outliers)

print(results, width = Inf)

# ── Save ──────────────────────────────────────────────────────────────────────
out_path <- file.path(GEMMA_DIR, "do_bidirectional_mr_results.csv")
write_csv(results, out_path)
cat("\nResults written to:", out_path, "\n")

# ── Instrument summary (save for supplement) ──────────────────────────────────
iv_summary <- bind_rows(
  iv_chol |> mutate(arm = "cholesterol → calcium",
                    exp_beta = beta_chol, exp_se = se_chol,
                    out_beta = beta_ca,   out_se  = se_ca),
  iv_ca   |> mutate(arm = "calcium → cholesterol",
                    exp_beta = beta_ca,   exp_se = se_ca,
                    out_beta = beta_chol, out_se  = se_chol)
) |>
  select(arm, snp, chr_num, pos, a1, a0, af,
         exp_beta, exp_se, out_beta, out_se,
         pval_chol, pval_ca, F_chol, F_ca)

write_csv(iv_summary,
          file.path(GEMMA_DIR, "do_bidirectional_mr_instruments.csv"))
cat("Instrument list written to: do_bidirectional_mr_instruments.csv\n")
