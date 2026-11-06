## Bidirectional MR: cholesterol ↔ calcium in DO mice (Svenson HFD v12)
##
## Same individuals provide both traits → 100% sample overlap.
##
## Co-primary methods:
##   MRBEE.IMRP.UV  — bias-corrected IVW; estimates error covariance (Rxy)
##                    from null SNPs in the GWAS data, which directly captures
##                    the overlap-induced correlation without needing phenotype files
##   IVW (overlap-corrected) — first-order bias correction to numerator and
##                    denominator; independent implementation for cross-check
##
## MR-APSS is NOT run: it requires LD scores for est_paras(), and no
## DO-mouse LD score reference panel is available.
##
## Run locally (where MRBEE is installed):
##   Rscript do_bidirectional_mr.R
##
## Install MRBEE if needed:
##   remotes::install_github("noahlorinczcomi/MRBEE")

suppressPackageStartupMessages({
  library(tidyverse)
  library(MRBEE)
})

# ── Paths (local copies of cluster files) ────────────────────────────────────
DATA_DIR <- here::here("Mouse Genetics/Cholesterol/GWAS")
# Fall back to script location if here::here() is not configured
if (!dir.exists(DATA_DIR)) {
  DATA_DIR <- normalizePath(
    file.path(dirname(sys.frame(1)$ofile), "."),
    mustWork = FALSE
  )
}

# ── Thresholds ───────────────────────────────────────────────────────────────
P_THRESH <- 5e-6   # instrument selection; relax to 1e-5 if < 3 instruments
WIND_MB  <- 1      # one top SNP per Mb window
MIN_F    <- 10     # minimum per-SNP F-statistic

# ── Load GWAS summary stats ───────────────────────────────────────────────────
read_gemma <- function(path, trait) {
  read_tsv(path, col_types = cols(.default = "d", rs = "c",
                                  allele1 = "c", allele0 = "c"),
           show_col_types = FALSE) |>
    rename(snp = rs, pos = ps, a1 = allele1, a0 = allele0, pval = p_wald) |>
    mutate(trait   = trait,
           chr_num = as.integer(str_extract(snp, "^[^_]+"))) |>
    filter(!is.na(beta), !is.na(se), se > 0)
}

chol <- read_gemma(file.path(DATA_DIR, "cholesterol_all.assoc.txt"), "cholesterol")
ca   <- read_gemma(file.path(DATA_DIR, "calcium2_all.assoc.txt"),    "calcium2")

cat("Cholesterol GWAS SNPs:", nrow(chol), "\n")
cat("Calcium GWAS SNPs:",     nrow(ca),   "\n")

# ── Harmonise: inner join on SNP id ──────────────────────────────────────────
harm <- inner_join(
  chol |> select(snp, chr_num, pos, a1, a0, af,
                 beta_chol = beta, se_chol = se, pval_chol = pval),
  ca   |> select(snp, beta_ca = beta, se_ca = se, pval_ca = pval),
  by = "snp"
) |>
  mutate(F_chol = (beta_chol / se_chol)^2,
         F_ca   = (beta_ca   / se_ca)^2)

cat("Harmonised SNPs:", nrow(harm), "\n")

# ── Rxy: MRBEE error covariance from null SNPs ────────────────────────────────
# Uses SNPs where both traits are non-significant to estimate the correlation
# induced by sharing the same individuals. This is MRBEE's recommended approach
# (Zhu et al. 2015 / MRBEE vignette) and avoids needing phenotype files.
Z_null <- harm |>
  filter(pval_chol > 0.05, pval_ca > 0.05) |>
  transmute(z_chol = beta_chol / se_chol,
            z_ca   = beta_ca   / se_ca) |>
  as.matrix()

Rxy <- tryCatch(
  as.matrix(MRBEE::errorCov(ZMatrix = Z_null)),
  error = function(e) {
    message("errorCov() failed (", e$message, "); falling back to cor(Z_null)")
    cor(Z_null, use = "complete.obs")
  }
)
cat("Rxy (error covariance):\n"); print(round(Rxy, 4))

# ── Instrument selection: top SNP per Mb window ──────────────────────────────
select_instruments <- function(df, pval_col, f_col,
                               p_thresh = P_THRESH,
                               wind_mb  = WIND_MB,
                               min_f    = MIN_F) {
  sig <- df |>
    filter(.data[[pval_col]] < p_thresh) |>
    arrange(.data[[pval_col]])

  if (nrow(sig) == 0) {
    warning("No SNPs below p < ", p_thresh, " for ", pval_col,
            ". Try relaxing P_THRESH to 1e-5.")
    return(sig[0, ])
  }

  kept <- list()
  while (nrow(sig) > 0) {
    top <- sig[1, ]
    kept[[length(kept) + 1]] <- top
    sig <- sig |>
      filter(!(chr_num == top$chr_num & abs(pos - top$pos) < wind_mb * 1e6))
  }
  out <- bind_rows(kept) |> filter(.data[[f_col]] >= min_f)
  cat(pval_col, "→", nrow(out), "instruments | mean F:",
      round(mean(out[[f_col]]), 1), "\n")
  out
}

iv_chol <- select_instruments(harm, "pval_chol", "F_chol")
iv_ca   <- select_instruments(harm, "pval_ca",   "F_ca")

# ── MRBEE.IMRP.UV ─────────────────────────────────────────────────────────────
run_mrbee <- function(ivs, exp_beta, exp_se, out_beta, out_se, Rxy, label) {
  if (nrow(ivs) < 3) {
    warning(label, ": fewer than 3 instruments, skipping MRBEE")
    return(NULL)
  }
  res <- tryCatch(
    MRBEE::MRBEE.IMRP.UV(
      by    = ivs[[out_beta]],
      bx    = ivs[[exp_beta]],
      byse  = ivs[[out_se]],
      bxse  = ivs[[exp_se]],
      Rxy   = Rxy
    ),
    error = function(e) { warning(label, " MRBEE failed: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)

  b  <- as.numeric(res$theta)[1]
  v  <- if (!is.null(res$vartheta)) res$vartheta else res$covtheta
  se <- sqrt(as.numeric(v)[1])
  n_pleio <- if (!is.null(res$delta)) sum(abs(as.numeric(res$delta)) > 1e-8) else NA_integer_

  tibble(arm = label, method = "MRBEE.IMRP.UV",
         n_ivs = nrow(ivs),
         mean_F = round(mean((ivs[[exp_beta]] / ivs[[exp_se]])^2), 1),
         beta = b, se = se,
         pval = 2 * pnorm(abs(b / se), lower.tail = FALSE),
         n_pleiotropic = n_pleio)
}

# ── Overlap-corrected IVW (independent cross-check) ──────────────────────────
# Corrects numerator (covariance bias) and denominator (winner's curse) using
# the off-diagonal of Rxy as ρ.
run_ivw_corrected <- function(ivs, exp_beta, exp_se, out_beta, out_se,
                              Rxy, label) {
  if (nrow(ivs) < 2) return(NULL)
  bx   <- ivs[[exp_beta]]; bxse <- ivs[[exp_se]]
  by   <- ivs[[out_beta]]; byse <- ivs[[out_se]]
  rho  <- Rxy[1, 2]
  w    <- 1 / byse^2

  num   <- sum(w * bx * by)   - rho * sum(w * bxse * byse)
  denom <- sum(w * bx^2)      - sum(w * bxse^2)

  if (denom <= 0) {
    warning(label, ": IVW corrected denominator ≤ 0 — instruments too weak")
    return(NULL)
  }
  beta <- num / denom
  se   <- sqrt(sum(w * (by - beta * bx)^2) / ((length(bx) - 1) * sum(w * bx^2)))

  tibble(arm = label, method = "IVW (overlap-corrected)",
         n_ivs = length(bx),
         mean_F = round(mean((bx / bxse)^2), 1),
         beta = beta, se = se,
         pval = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
         n_pleiotropic = NA_integer_)
}

# ── Run bidirectional MR ──────────────────────────────────────────────────────
cat("\n--- Cholesterol → Calcium ---\n")
r1 <- run_mrbee(iv_chol, "beta_chol", "se_chol", "beta_ca", "se_ca",
                Rxy, "cholesterol → calcium")
r2 <- run_ivw_corrected(iv_chol, "beta_chol", "se_chol", "beta_ca", "se_ca",
                        Rxy, "cholesterol → calcium")

cat("\n--- Calcium → Cholesterol ---\n")
r3 <- run_mrbee(iv_ca, "beta_ca", "se_ca", "beta_chol", "se_chol",
                Rxy, "calcium → cholesterol")
r4 <- run_ivw_corrected(iv_ca, "beta_ca", "se_ca", "beta_chol", "se_chol",
                        Rxy, "calcium → cholesterol")

results <- bind_rows(r1, r2, r3, r4) |>
  mutate(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se) |>
  relocate(arm, method, n_ivs, mean_F, beta, se, ci_lo, ci_hi, pval, n_pleiotropic)

cat("\n=== Summary ===\n")
print(results, width = Inf)

# ── Save ──────────────────────────────────────────────────────────────────────
write_csv(results,
          file.path(DATA_DIR, "do_bidirectional_mr_results.csv"))

write_csv(
  bind_rows(
    iv_chol |> mutate(arm = "cholesterol → calcium"),
    iv_ca   |> mutate(arm = "calcium → cholesterol")
  ),
  file.path(DATA_DIR, "do_bidirectional_mr_instruments.csv")
)

cat("\nResults       → do_bidirectional_mr_results.csv\n")
cat("Instruments   → do_bidirectional_mr_instruments.csv\n")
cat("Rxy off-diag  →", round(Rxy[1, 2], 4),
    "(overlap correlation; should match phenotypic r)\n")

# ── Combine with MR-APSS results if available ────────────────────────────────
apss_path <- file.path(DATA_DIR, "do_mr_apss_results.csv")
if (file.exists(apss_path)) {
  apss <- read_csv(apss_path, show_col_types = FALSE)
  combined <- bind_rows(results, apss) |>
    arrange(arm, method)
  write_csv(combined, file.path(DATA_DIR, "do_bidirectional_mr_all_methods.csv"))
  cat("\n=== Combined (MRBEE + MR-APSS) ===\n")
  print(combined |> select(arm, method, n_ivs, beta, se, pval), width = Inf)
  cat("\nCombined → do_bidirectional_mr_all_methods.csv\n")
} else {
  cat("\nNote: do_mr_apss_results.csv not found. Run do_mr_apss.R on the cluster,\n")
  cat("download the result to this folder, then re-run this script to combine.\n")
}
