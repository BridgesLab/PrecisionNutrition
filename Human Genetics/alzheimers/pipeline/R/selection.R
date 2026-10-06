# selection.R — survival-collider (SlopeHunter) inputs for the pipeline.
#
# b_SH is NOT refit by default: the pipeline applies the fits Step 2 (exposure_ad_indexevent.qmd)
# already uses, with the same `trusted` gate (Kunkle's fit is not identified, so it is never
# applied). What this module adds, now that the full Bellenguez sumstats are local, is the check
# Step 2 left open (its caveat 0): is Bellenguez's b_SH = -0.944 distinguishable from what
# SlopeHunter returns on permuted (signal-free) data? Engine: R/fit_selection_slope.R.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# Step 2's b_SH table, same files and gate. One row per AD dataset key.
load_bsh_table <- function(bsh_cfg) {
  purrr::imap(bsh_cfg, function(x, key) {
    fit <- if (file.exists(x$fits_csv)) readr::read_csv(x$fits_csv, show_col_types = FALSE) |>
      filter(set == x$fit_set) |> slice(1) else NULL
    tibble(dataset = key, b_SH = fit$b_SH %||% NA_real_, se_SH = fit$se_SH %||% NA_real_,
           trusted = isTRUE(x$trusted) && !is.null(fit), reason = x$reason,
           source = x$fits_csv)
  }) |> bind_rows()
}

# AD x lifespan merge in the layout fit_selection_slope.R's diagnostics expect, built from the
# pipeline's standard-format tables (AD aligned to the lifespan effect allele), then LD-clumped on
# the lifespan p-value exactly as the Kunkle fit was (p1 1e-3, p2 0.01, r2 0.01, 1 Mb).
selection_merge <- function(g_ad, g_life, bfile, plink2, life_p_keep = 0.01) {
  life <- g_life[p < life_p_keep]
  m <- merge(life[, .(SNP, chr, pos, ea_l = ea, oa_l = oa, beta_life = beta, se_life = se, p_life = p)],
             g_ad[, .(SNP, ea_a = ea, oa_a = oa, beta_ad = beta, se_ad = se, p_ad = p)], by = "SNP")
  m <- m[!paste0(ea_l, oa_l) %chin% c("AT", "TA", "CG", "GC")]
  same <- m$ea_a == m$ea_l & m$oa_a == m$oa_l
  flip <- m$ea_a == m$oa_l & m$oa_a == m$ea_l
  m <- m[same | flip]
  m[, beta_ad_aligned := fifelse(ea_a == ea_l, beta_ad, -beta_ad)]
  m[, apoe := chr == "19" & pos > 43.9e6 & pos < 46.5e6]   # spans both builds' APOE windows
  # Clump on lifespan (p1 1e-3 keeps the SNPs SlopeHunter fits on; p2 0.01 assigns the rest).
  f <- tempfile("selclump_")
  fwrite(m[, .(ID = SNP, P = p_life)], paste0(f, ".tsv"), sep = "\t")
  status <- system2(plink2, c("--bfile", shQuote(bfile), "--clump", shQuote(paste0(f, ".tsv")),
                              "--clump-id-field", "ID", "--clump-p-field", "P", "--clump-p1", "1e-3",
                              "--clump-p2", "0.01", "--clump-r2", "0.01", "--clump-kb", "1000",
                              "--out", shQuote(f), "--silent"), stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("plink2 clump for the selection merge failed")
  idx <- fread(paste0(f, ".clumps"), select = "ID")$ID
  m[SNP %chin% idx]
}

# Identifiability, the permutation null and a refit, on the APOE-excluded clumped set.
bsh_checks <- function(m, applied_b_sh, n_perm = 8) {
  prim <- m[apoe == FALSE]
  ident <- bsh_identifiability(prim)
  null <- bsh_null_permutation(prim, n_perm = n_perm)
  refit <- fit_slopehunter_one(prim, "excl_APOE (pipeline refit)")
  tibble(n_fit = ident$n_fit, pearson_r = ident$pearson_r, sign_agreement = ident$sign_agreement,
         identifiability = ident$verdict,
         b_SH_applied = applied_b_sh, b_SH_refit = refit$b_SH, se_SH_refit = refit$se_SH,
         null_median = null$null_median, null_min = null$null_min, null_max = null$null_max,
         null_sd = null$null_sd,
         # Observed is distinguishable from the no-signal output only if it sits outside the
         # permutation range by more than the bootstrap SE.
         distinguishable_from_null = (applied_b_sh < null$null_min - 2 * refit$se_SH) |
                                     (applied_b_sh > null$null_max + 2 * refit$se_SH))
}
