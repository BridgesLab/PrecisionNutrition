# pipeline.R — turn config.yml into the dataset and association tables _targets.R maps over, and
# assemble the final tables.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr) })

dataset_table <- function(cfg) {
  imap(cfg$datasets, \(d, key) tibble(
    key = key, id = d$id, label = d$label, source = d$source, file = d$file %||% NA_character_,
    n = d$sample_size, type = d$type, cohorts = list(d$cohorts),
    adjusted_for = d$adjusted_for %||% NA_character_, p_from_beta = isTRUE(d$p_from_beta))) |>
    bind_rows()
}

# One row per association x outcome arm, with the pre-specified methods (config `analysis:`):
#   main (MR-APSS) and sensitivity (CAUSE) for every association; none for classical_only blocks.
# `overlap` (shared cohorts, incl. the lifespan GWAS's for the SlopeHunter arm, whose adjusted
# outcome is built from it) is kept for information; MR-APSS and CAUSE both model overlap.
# Until 2026-10-05 the rule was MRAID without overlap / CAUSE with it (see config.yml).
association_table <- function(cfg, ds) {
  coh <- setNames(ds$cohorts, ds$key)
  trusted <- names(keep(cfg$slopehunter$fits, \(f) isTRUE(f$trusted)))
  life <- cfg$slopehunter$selection_axis
  fwd <- imap(cfg$associations, \(blk, block) expand_grid(exposure = blk$exposures, outcome = blk$outcomes) |>
                mutate(block = block, classical_only = isTRUE(blk$classical_only))) |>
    bind_rows() |>
    filter(exposure != outcome)
  # Reverse directions (analysis$reverse_directions): every non-classical-only pair, swapped. The
  # reversed outcome is never a SlopeHunter-trusted AD GWAS, so these get the raw arm only.
  rev <- if (isTRUE(cfg$analysis$reverse_directions))
    fwd |> filter(!classical_only) |> distinct(exposure, outcome, block) |>
      # Swap through a temporary column: transmute() evaluates in order, so a direct swap would
      # read the already-overwritten exposure.
      transmute(tmp = exposure, exposure = outcome, outcome = tmp, block = paste0(block, "_reverse"),
                classical_only = FALSE) |>
      select(-tmp) |>
      anti_join(fwd, by = c("exposure", "outcome"))
  bind_rows(fwd, rev) |>
    mutate(arm = map(outcome, \(o) c("raw", if (o %in% trusted) "slopehunter"))) |>
    unnest(arm) |>
    mutate(
      outcome_data = if_else(arm == "slopehunter", paste0(outcome, "_sh"), outcome),
      overlap = pmap_chr(list(exposure, outcome, arm), \(e, o, a) {
        shared <- intersect(coh[[e]], c(coh[[o]], if (a == "slopehunter") coh[[life]]))
        if (length(shared)) paste(shared, collapse = "+") else "none"
      }),
      primary_method = if_else(classical_only, "none", cfg$analysis$main),
      sensitivity_method = if_else(classical_only, "none", cfg$analysis$sensitivity),
      assoc_id = paste(exposure, outcome, arm, sep = "__"))
}

# Final results table: one row per association x method, classical suite and primary method
# together; diagnostics (Q, I2, I2_GX, R2, F, Steiger) repeat on every classical row of an
# association. `is_primary` marks the pre-specified primary estimate.
assemble_results <- function(assoc, classical_rows, primary_rows, stability) {
  bind_rows(classical_rows, primary_rows) |>
    left_join(assoc |> select(assoc_id, block, primary_method, sensitivity_method, overlap), by = "assoc_id") |>
    mutate(is_primary = method == primary_method,
           is_sensitivity = method == sensitivity_method,
           ci_lo = if ("ci_lo" %in% names(pick(everything()))) coalesce(ci_lo, b - 1.96 * se) else b - 1.96 * se,
           ci_hi = if ("ci_hi" %in% names(pick(everything()))) coalesce(ci_hi, b + 1.96 * se) else b + 1.96 * se) |>
    left_join(stability, by = c("assoc_id", "method")) |>
    # Overlap, to filter overlap-naive methods: the configured shared cohorts (`overlap`), MR-APSS's
    # empirical cross-trait LDSC intercept per pair, and whether each method models overlap.
    left_join(primary_rows |> filter(method == "MR-APSS") |> select(assoc_id, overlap_C12 = any_of("apss_C12")),
              by = "assoc_id") |>
    mutate(method_handles_overlap = method %in% c("MR-APSS", "MR-APSS (p<5e-8)", "CAUSE", "MRBEE")) |>
    relocate(assoc_id, exposure, outcome, arm, method, is_primary, is_sensitivity, scale, b, se,
             ci_lo, ci_hi, p, n_snp) |>
    arrange(block, exposure, outcome, arm, desc(is_primary), desc(is_sensitivity), method)
}

write_csv_target <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(x, path)
  path
}

# ---- per-association wrappers (what the targets call) ------------------------------------------
# Each returns list(rows, instruments, attrition) so the final tables are plain bind_rows().

meta_cols <- function(assoc_id, exposure, outcome, arm)
  tibble(assoc_id = assoc_id, exposure = exposure, outcome = outcome, arm = arm)

analyse_classical <- function(inst_sel, g_exp, g_out, assoc_id, exposure, outcome, arm,
                              presso_nboot, seed) {
  h <- harmonise_pair(inst_sel$inst, g_out, assoc_id, "classical")
  att <- bind_rows(inst_sel$attrition |> mutate(association = assoc_id), h$attrition)
  if (is.null(h$dat) || nrow(h$dat) < 3)
    return(list(rows = bind_cols(meta_cols(assoc_id, exposure, outcome, arm),
                                 tibble(method = "all", note = "not run: < 3 harmonised SNPs")),
                instruments = NULL, attrition = att))
  Rxy <- error_cor(g_exp, g_out)
  diag <- bind_cols(instrument_diagnostics(h$dat) |> rename(n_snp_diag = n_snp), steiger_z(h$dat))
  rows <- run_classical_suite(h$dat, Rxy = Rxy, presso_nboot = presso_nboot, seed = seed)
  list(rows = bind_cols(meta_cols(assoc_id, exposure, outcome, arm), rows, diag),
       instruments = as_tibble(h$dat) |>
         transmute(assoc_id = assoc_id, instrument_set = "classical", SNP,
                   effect_allele = effect_allele.exposure, other_allele = other_allele.exposure,
                   beta_exposure = beta.exposure, se_exposure = se.exposure,
                   p_exposure = pval.exposure, beta_outcome = beta.outcome,
                   se_outcome = se.outcome, p_outcome = pval.outcome),
       attrition = att)
}

analyse_mraid <- function(inst_sel, g_out, assoc_id, exposure, outcome, arm, bim, bfile, plink2,
                          n_exp, n_out, cfg, max_candidates = Inf) {
  # Pre-specified rule (config analysis$mraid_rule): MRAID was unreliable with thousands of
  # correlated candidates (lean mass, fat-free mass), so above the cap it is recorded as not run.
  n_cand <- nrow(inst_sel$inst)
  if (n_cand > max_candidates)
    return(list(rows = bind_cols(meta_cols(assoc_id, exposure, outcome, arm),
                                 tibble(method = "MRAID", scale = "SD(outcome) per SD(exposure)",
                                        n_snp = as.integer(n_cand),
                                        note = sprintf("not run: %d candidates > %d (pre-specified cap)",
                                                       n_cand, as.integer(max_candidates)))),
                instruments = NULL,
                attrition = inst_sel$attrition |> mutate(association = assoc_id)))
  r <- run_mraid(inst_sel$inst, g_out, assoc_id, bim, bfile, plink2, n_exp, n_out, cfg)
  list(rows = bind_cols(meta_cols(assoc_id, exposure, outcome, arm), r$result,
                        tibble(mraid_r2_used = inst_sel$r2_used)),
       instruments = if (!is.null(r$instruments))
         as_tibble(r$instruments) |> mutate(assoc_id = assoc_id, instrument_set = "mraid", .before = 1),
       attrition = bind_rows(inst_sel$attrition |> mutate(association = assoc_id), r$attrition))
}

analyse_cause <- function(g_exp, g_out, assoc_id, exposure, outcome, arm, hm3, n_exp, n_out,
                          bfile, plink2, cfg, seed) {
  r <- run_cause(g_exp, g_out, assoc_id, hm3, n_exp, n_out, bfile, plink2, cfg, seed)
  list(rows = bind_cols(meta_cols(assoc_id, exposure, outcome, arm), r$result),
       instruments = as_tibble(r$instruments) |> mutate(assoc_id = assoc_id, instrument_set = "cause", .before = 1),
       attrition = r$attrition, weights = r$weights |> mutate(assoc_id = assoc_id, .before = 1),
       elpd = r$elpd |> mutate(assoc_id = assoc_id, .before = 1),
       posterior = r$posterior)
}

# Step 2's missing term, filled: E -> lifespan (IVW-MRE, classical instruments) and the
# estimate-level correction naive - b_SH * (E -> lifespan), exactly Step 2's formula, for every
# exposure into every trusted-fit AD outcome. Compare with the SNP-level "slopehunter" arm rows.
step2_gap_fill <- function(results, bsh_table) {
  ivw <- results |> filter(method == "IVW-MRE")
  life <- ivw |> filter(outcome == "lifespan") |> select(exposure, b_life = b, se_life = se, p_life = p)
  d <- ivw |> filter(arm == "raw", outcome %in% bsh_table$dataset[bsh_table$trusted]) |>
    select(exposure, outcome, b_naive = b, se_naive = se, p_naive = p) |>
    left_join(life, by = "exposure") |>
    left_join(bsh_table |> select(outcome = dataset, b_SH, se_SH), by = "outcome")
  bind_cols(d, correct_estimate(d$b_naive, d$se_naive, d$b_life, d$se_life, d$b_SH, d$se_SH) |>
              select(b_corr, se_corr, ci_lo, ci_hi, p_corr))
}
