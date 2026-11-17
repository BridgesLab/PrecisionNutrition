# Where the files in `results/pipeline/` came from

**Stamps (from 2026-11-17).** Every CSV written by the pipeline (`write_csv_target()`) or by
`pipeline/bsh_propagate_existing.R` now gets a stamp file at `stamps/<file>.txt` beside it
(`pipeline/R/provenance.R`). The stamp records the file's md5, the step that wrote it, the time and
host, the git commit (or "not a git checkout" on a Great Lakes copy), a digest of the pipeline code,
and package versions. A file with no stamp predates this. Its history is the table below,
reconstructed from `git log` on 2026-11-17. **Add a row to the table whenever a run or
post-processing step changes a file here.**

Every CSV is written by a `targets` target through `write_csv_target()`. The target names are below.
Per `pipeline/bsh_propagate_existing.R`, the CSVs in git come from Great Lakes runs
(`pipeline/run_greatlakes.sbatch`). The local `_targets/` store is older, so don't rebuild these
files locally with a bare `tar_make()`.

## History

| Commit | Date | What changed | Files |
|---|---|---|---|
| `e3f7e59` | 2026-10-06 | First pipeline results | `results.csv`, `instruments_{primary,classical}.csv`, `instrument_attrition.csv.gz`, `associations.csv`, `datasets.csv`, `mrapss_background.csv`, `step2_lifespan_gap_fill.csv`, `bsh_null_checks.csv` |
| `42568d8` | 2026-10-08 | **Run of 2026-10-07**, pipeline code at `071fa52` (adds latent confounding, MR-APSS p < 5×10⁻⁸ sensitivity, reverse directions, fat mass). Source of every univariable number in `SUMMARY.md` | `results.csv`, `instruments_*.csv`, `instrument_attrition.csv.gz`, `associations.csv`, `datasets.csv`, `latent_confounding.csv`, `mrapss_background.csv`, `step2_lifespan_gap_fill.csv` |
| `38987ad` | 2026-10-12 | Fasting vs 2 h glucose contrast, **Step 0 (outcome-blind)**, plus the frozen `assignments/` | `contrast/step0_{strength,overlap,instruments,attrition,cis_instruments,cis_summary}.csv` |
| `cccee6b` | 2026-10-13 | Contrast outcome steps (Bellenguez raw + SlopeHunter, Kunkle) | `contrast/{results,cluster_results,drug_target_results,error_correlations,instrument_outcome_data}.csv` |
| `022604a` | 2026-10-15 | Drug-target follow-up; post-freeze class (bromocriptine) | `contrast/drug_target_snps.csv`, `contrast/post_freeze_*.csv` |
| `6814fd5` | 2026-10-29 | Uffelmann 2026 new-cohorts and clinical-only outcomes added (post-freeze); LD map of Xue instruments onto Suzuki clusters | `contrast/{results,cluster_results,drug_target_results,drug_target_snps,error_correlations,instrument_outcome_data,xue_cluster_map}.csv` |
| `6930e3e` | 2026-10-30 | Combined Uffelmann (all cohorts) sensitivity outcome | same `contrast/` files as `6814fd5`, except `xue_cluster_map.csv` |
| `7d9cad5` | 2026-11-11 | **Post-processing, not a pipeline run:** `pipeline/bsh_propagate_existing.R` added the b_SH uncertainty to the SlopeHunter rows in place (`se`, `p`, CI; conditional values kept in `*_cond_bsh`; new `p_bsh_cons`, `bsh_tip`) and wrote `bsh_grid.csv` | `results.csv`, `bsh_grid.csv`, `contrast/{results,cluster_results,drug_target_results,post_freeze_drug_target_results}.csv` |
| (2026-11-17 commit) | 2026-11-17 | Local, by name: `tar_make(names = c("datasets_used", "csv_datasets"))`. The file had last been written on 2026-10-08, before `fg_chen`, `t2d_suzuki_eur` and the three Uffelmann files were added to the config. Existing rows unchanged; five rows added. First file with a stamp | `datasets.csv`, `stamps/datasets.csv.txt` |

`figures/<assoc_id>/` (scatter, funnel and leave-one-out plots per association) came from the same
runs as `results.csv` (targets `diag_classical`, `diag_apss`).

## Which target writes each file

| File | Target |
|---|---|
| `results.csv` | `csv_results` (with `propagate_bsh()`) |
| `instruments_primary.csv`, `instruments_classical.csv` | `csv_instruments`, `csv_instruments_classical` |
| `instrument_attrition.csv.gz` | `csv_attrition` |
| `associations.csv`, `datasets.csv` | `csv_associations`, `csv_datasets` |
| `mrapss_background.csv`, `latent_confounding.csv` | `csv_apss_background`, `csv_latent` |
| `step2_lifespan_gap_fill.csv`, `bsh_null_checks.csv`, `bsh_grid.csv` | `csv_step2_gap`, `csv_bsh_checks`, `csv_bsh_grid` |
| `contrast/step0_*.csv` | `csv_contrast_{strength,overlap,instruments,attrition}`, `csv_cis_{summary,instruments}` |
| `contrast/results.csv`, `cluster_results.csv` | `csv_contrast_results`, `csv_contrast_groups` |
| `contrast/instrument_outcome_data.csv`, `error_correlations.csv` | `csv_contrast_data`, `csv_contrast_errcor` |
| `contrast/drug_target_results.csv`, `drug_target_snps.csv` | `csv_cis_results`, `csv_cis_snps` |
| `contrast/xue_cluster_map.csv` | `csv_xue_cluster_map` |
| `contrast/post_freeze_*.csv` | `csv_pf` |

## Inputs from outside this folder

| Input | Where | Used for |
|---|---|---|
| Bellenguez b_SH fit | `R/fits/bellenguez_lifespan_slopehunter_fits.csv` (set `excl_APOE (PRIMARY)`), from `scripts/refit_bellenguez_bsh.R` since 2026-11-17. Runs up to and including `7d9cad5` read `smoking/results/indexevent_slopehunter_fits.csv`, whose PRIMARY row is identical to full precision | every SlopeHunter-corrected row |
| Kunkle and Uffelmann b_SH fits | `R/fits/*_lifespan_slopehunter_fits.csv` (marked `trusted: false`, never applied) | the identification check only |
| HapMap3 SNP list | `../calcium-cholesterol/raw_data/reference/w_hm3.snplist` | MRAID and CAUSE instrument sets |
| LD scores | `../calcium-cholesterol/raw_data/reference/eur_w_ld_chr` | MR-APSS background |
| 1000G EUR panel | `data/cache/EUR.{bed,bim,fam}` (gitignored) | clumping, LD |
| GWAS summary statistics | `data/cache/` (gitignored), sources in `pipeline/config.yml` → `datasets:` | everything |
| Suzuki 2024 clusters | `data/cache/suzuki2024/suzuki2024_supplementary_tables.xlsx`, sheet ST6 | frozen `assignments/` |
| Open Targets T2D drug targets | frozen snapshot `assignments/opentargets_t2d_drug_targets.csv` | drug-target classes |

## Going forward

- Stamps are written automatically. Add one row per run or post-processing step to the history table.
- Changing `write_csv_target()` (2026-11-17) makes every `csv_*` target outdated. The next full
  Great Lakes run rewrites and stamps them all. The estimates do not change, because the upstream
  targets are untouched.
- Locally, keep building by name. A bare `tar_make()` would rebuild from the older local store.
