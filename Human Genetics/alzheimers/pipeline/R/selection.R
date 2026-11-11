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

# ---- b_SH uncertainty: propagation and sensitivity grid -----------------------------------------
# adjust_outcome_gw() treats b_SH as known: the adjusted SE adds only b_SH^2 * se_life^2. The error
# in b_SH is shared by every SNP, so it cannot be added SNP by SNP (IVW would treat it as independent
# noise and shrink it by ~1/J). Zhang et al. 2026 (AJHG; medRxiv 2026.02.26.26347193) show that
# leaving it out gives SlopeHunter ~50% Type I error under the null. It is propagated per estimate:
#
#   theta(b) = theta_raw + b * d,   d = dtheta/db = (theta_SH - theta_raw) / b_SH
#   se_total^2 = se_cond^2 + d^2 * se_SH^2
#
# d is read off the raw (b = 0) and adjusted (b = b_SH) rows of the same estimator. For IVW it is
# minus the estimate on the lifespan GWAS (up to weights); for MR-APSS, CAUSE, MRAID, MRBEE, PRESSO
# and the median/mode it is a secant approximation. se_cond is the SE given b_SH, as the methods
# report it. CIs widen on each side in quadrature (keeps CAUSE's asymmetry) and CAUSE's
# P(gamma < 0) is rescaled on the normal-quantile scale by se_cond / se_total.
#
# se_SH is SlopeHunter's bootstrap SE. Bellenguez's permutation-null SD (bsh_checks) is ~3x larger,
# so p_bsh_cons repeats the calculation with the larger of the two. bsh_tip is the b at which the
# estimate crosses 0 (compare with b_SH and the permutation-null median).

# Row keys (besides outcome and arm) that pair each adjusted row with its raw row, per table.
BSH_KEYS <- list(
  results     = c("exposure", "method"),
  contrast    = c("set", "step", "estimator", "term", "scale"),
  cluster     = c("set", "model", "group"),
  drug_target = c("set", "scale", "level", "target", "trait"))

BSH_ADDED <- c("se_cond_bsh", "p_cond_bsh", "ci_lo_cond_bsh", "ci_hi_cond_bsh",
               "p_gamma_below_0_cond_bsh", "dtheta_dbsh", "se_bsh_term", "p_bsh_cons", "bsh_tip")

# One row per trusted outcome: the applied b_SH, its bootstrap SE and the permutation null.
bsh_params <- function(bsh_table, bsh_checks) {
  bsh_table |> filter(trusted) |>
    select(outcome = dataset, b_SH, se_SH) |>
    left_join(bsh_checks |> select(outcome = dataset, null_median, null_sd), by = "outcome") |>
    mutate(se_SH_cons = pmax(se_SH, null_sd, na.rm = TRUE))
}

# Undo an earlier propagate_bsh(), so re-applying it never compounds the extra variance.
unpropagate_bsh <- function(tbl) {
  if (!"dtheta_dbsh" %in% names(tbl)) return(tbl)
  back <- c(se = "se_cond_bsh", p = "p_cond_bsh", ci_lo = "ci_lo_cond_bsh", ci_hi = "ci_hi_cond_bsh",
            p_gamma_below_0 = "p_gamma_below_0_cond_bsh")
  done <- !is.na(tbl$dtheta_dbsh)
  for (col in names(back)) if (all(c(col, back[[col]]) %in% names(tbl)))
    tbl[[col]][done] <- tbl[[back[[col]]]][done]
  if ("ci_excludes_0" %in% names(tbl)) {
    j <- done & !is.na(tbl$ci_excludes_0)
    tbl$ci_excludes_0[j] <- tbl$ci_lo[j] > 0 | tbl$ci_hi[j] < 0
  }
  tbl |> select(-any_of(BSH_ADDED))
}

# SlopeHunter rows joined to their raw rows and the outcome's b_SH. Internal.
bsh_pairs <- function(tbl, keys, par) {
  by <- c(keys, "outcome")
  raw <- tbl |> filter(arm == "raw") |> select(all_of(by), .b_raw = b, .se_raw = se)
  if (anyDuplicated(raw[by])) stop("bsh: raw rows are not unique on ", paste(by, collapse = ", "))
  tbl |> mutate(.row = row_number()) |> filter(arm == "slopehunter") |>
    left_join(raw, by = by) |>
    left_join(par |> select(outcome, .b_SH = b_SH, .se_SH = se_SH, .se_cons = se_SH_cons,
                            .null_median = null_median), by = "outcome") |>
    mutate(.slope = (b - .b_raw) / .b_SH)
}

# Replace se / p / CI of the SlopeHunter rows with values that include the b_SH uncertainty; keep
# the conditional values in *_cond_bsh. Raw rows are untouched. Idempotent.
propagate_bsh <- function(tbl, keys, par) {
  tbl <- unpropagate_bsh(tbl)
  if (!nrow(tbl) || !any(tbl$arm == "slopehunter")) return(tbl)
  has_ci <- all(c("ci_lo", "ci_hi") %in% names(tbl))
  has_pg <- "p_gamma_below_0" %in% names(tbl)
  z <- stats::qnorm(0.975)
  s <- bsh_pairs(tbl, keys, par) |> filter(is.finite(.slope))
  se_eq <- if (has_ci) dplyr::coalesce(s$se, (s$ci_hi - s$ci_lo) / (2 * z)) else s$se
  tau   <- abs(s$.slope) * s$.se_SH
  se_t  <- sqrt(se_eq^2 + tau^2)
  se_tc <- sqrt(se_eq^2 + (abs(s$.slope) * s$.se_cons)^2)
  i <- s$.row
  added <- setdiff(BSH_ADDED, c(if (!has_ci) c("ci_lo_cond_bsh", "ci_hi_cond_bsh"),
                                if (!has_pg) "p_gamma_below_0_cond_bsh"))
  for (col in added) tbl[[col]] <- NA_real_
  tbl$se_cond_bsh[i] <- s$se
  tbl$p_cond_bsh[i]  <- s$p
  tbl$dtheta_dbsh[i] <- s$.slope
  tbl$se_bsh_term[i] <- tau
  tbl$p_bsh_cons[i]  <- ifelse(is.na(s$p), NA_real_, 2 * stats::pnorm(-abs(s$b / se_tc)))
  tbl$bsh_tip[i]     <- ifelse(s$.slope == 0, NA_real_, -s$.b_raw / s$.slope)
  tbl$se[i] <- ifelse(is.na(s$se), NA_real_, se_t)
  tbl$p[i]  <- ifelse(is.na(s$p), NA_real_, 2 * stats::pnorm(-abs(s$b / se_t)))
  if (has_ci) {
    tbl$ci_lo_cond_bsh[i] <- s$ci_lo
    tbl$ci_hi_cond_bsh[i] <- s$ci_hi
    tbl$ci_lo[i] <- s$b - sqrt((s$b - s$ci_lo)^2 + (z * tau)^2)
    tbl$ci_hi[i] <- s$b + sqrt((s$ci_hi - s$b)^2 + (z * tau)^2)
    if ("ci_excludes_0" %in% names(tbl)) {     # CAUSE's decision rule; NA for the other methods
      j <- i[!is.na(tbl$ci_excludes_0[i])]
      tbl$ci_excludes_0[j] <- tbl$ci_lo[j] > 0 | tbl$ci_hi[j] < 0
    }
  }
  if (has_pg) {
    pg <- pmin(pmax(s$p_gamma_below_0, 1e-12), 1 - 1e-12)
    tbl$p_gamma_below_0_cond_bsh[i] <- s$p_gamma_below_0
    tbl$p_gamma_below_0[i] <- ifelse(is.na(pg), NA_real_, stats::pnorm(stats::qnorm(pg) * se_eq / se_t))
  }
  tbl
}

# Sensitivity grid: each SlopeHunter estimate at fixed values of b, interpolated between its raw
# (b = 0) and adjusted (b = b_SH) rows. SEs are conditional on that b, interpolated the same way.
# Rows without an SE (CAUSE) are left out.
bsh_sensitivity_grid <- function(tbl, keys, par, table) {
  tbl <- unpropagate_bsh(tbl)
  if (!nrow(tbl) || !any(tbl$arm == "slopehunter")) return(tibble())
  pts <- par |>
    transmute(outcome, `b = 0 (uncorrected)` = 0, `permutation-null median` = null_median,
              `b_SH - 1.96 SE` = b_SH - 1.96 * se_SH, `b_SH (applied)` = b_SH,
              `b_SH + 1.96 SE` = b_SH + 1.96 * se_SH) |>
    tidyr::pivot_longer(-outcome, names_to = "b_point", values_to = "b_value")
  s <- bsh_pairs(tbl, keys, par) |> filter(is.finite(.slope), !is.na(se), !is.na(.se_raw))
  s$row_key <- do.call(paste, c(unname(as.list(s[keys])), sep = " | "))
  s |> select(outcome, row_key, b_sh = b, se_sh = se, .b_raw, .se_raw, .b_SH) |>
    inner_join(pts, by = "outcome", relationship = "many-to-many") |>
    mutate(frac = b_value / .b_SH,
           estimate = .b_raw + frac * (b_sh - .b_raw),
           se = pmax(.se_raw + frac * (se_sh - .se_raw), pmin(.se_raw, se_sh)),
           p = 2 * stats::pnorm(-abs(estimate / se)),
           across(c(b_value, estimate, se, p), \(x) signif(x, 10)),
           table = table) |>
    select(table, outcome, row_key, b_point, b_value, estimate, se, p)
}
