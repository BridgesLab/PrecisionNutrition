# Figure and Table Key

This document shows where each figure panel and table comes from: the script, the chunk that draws
it, the output file, and the data that chunk reads. All code and data are in
**`Human Genetics/alzheimers/`**, and every path below is relative to that folder. Figures are made
in quarto and assembled in Adobe Illustrator.

**Numbering lives only here.** Output files are named by content (`figures/ms-<chunk label>-1.pdf`),
so renumbering, or reformatting for a journal, changes this table and nothing else. The numbering
below is a **draft** for the bioRxiv version.

**Workflow.** Dave writes `Manuscript.docx`; Claude never edits it. After each section is
drafted, Claude checks every method statement and number in it against the code and CSVs, traces
each back to its source (the table below, then `results/pipeline/RUN_INFO.md` and the stamps), and
replies with discrepancies and suggested wording. The STROBE-MR checklist and this table are updated
to match. Figures are assembled in Illustrator as the text is written. Place the `figures/ms-*.pdf`
files as **linked**, not embedded, so a re-render updates the artwork, and keep the Figure column
below in step with the Illustrator layout.

**Publication vs talk figures.** `manuscript_figures.qmd` makes the greyscale print figures (`ms-`).
`t2d_ad_figures.qmd` makes colour, large-font versions of the same panels for talks (`talk-`). Both
read the same CSVs, so the numbers agree.

## Main figures

| Figure | Content | Script (chunk) | File | Reads |
| ------ | ------- | -------------- | ---- | ----- |
| Figure 1 | Study design: IV assumptions and the survival collider | `mr_dags.qmd` (`overall`, `collider`) | `figures/dags-overall-1.pdf`, `figures/dags-collider-1.pdf` | nothing (illustrative) |
| Figure 2A | Exposure and AD-liability effects on parental lifespan (both arms of the collider) | `manuscript_figures.qmd` (`survival-collider`) | `figures/ms-survival-collider-1.pdf` | `results/pipeline/results.csv`; `ad_longevity_arm.csv` (`exposure_survival_screen.qmd`) |
| Figure 2B | T2D, BMI, lean mass → AD, uncorrected vs SlopeHunter-corrected | `manuscript_figures.qmd` (`ad-body-correction`) | `figures/ms-ad-body-correction-1.pdf` | `results/pipeline/results.csv` |
| Figure 3A | T2D, fasting and 2 h glucose → AD, both arms | `manuscript_figures.qmd` (`ad-glycemic`) | `figures/ms-ad-glycemic-1.pdf` | `results/pipeline/results.csv` |
| Figure 3B | Fasting vs 2 h glucose: MVMR and the pre-specified contrast regression, Suzuki and Xue instruments | `manuscript_figures.qmd` (`fg-2h`) | `figures/ms-fg-2h-1.pdf` | `results/pipeline/contrast/results.csv` |
| Figure 3C | Fasting vs 2 h glucose in Kunkle and Uffelmann 2026 | `manuscript_figures.qmd` (`fg-2h-replication`) | `figures/ms-fg-2h-replication-1.pdf` | `results/pipeline/contrast/results.csv` |
| Figure 4 | T2D liability → AD by Suzuki 2024 cluster, beside each cluster's glucose effects | `manuscript_figures.qmd` (`clusters`) | `figures/ms-clusters-1.pdf` | `contrast/cluster_results.csv`; glucose slopes computed in the chunk from `contrast/instrument_outcome_data.csv` |
| Figure 5A | Glucose-lowering drug-target classes → AD, beside their glucose effects | `manuscript_figures.qmd` (`drug-targets`) | `figures/ms-drug-targets-1.pdf` | `contrast/drug_target_results.csv` |
| Figure 5B | From T2D liability to the KATP channel (post hoc) | `manuscript_figures.qmd` (`katp-specificity`) | `figures/ms-katp-specificity-1.pdf` | computed in the chunk from `contrast/instrument_outcome_data.csv`, `contrast/xue_cluster_map.csv`, `contrast/step0_strength.csv` |

## Supplementary figures

| Figure | Content | Script (chunk) | File | Reads |
| ------ | ------- | -------------- | ---- | ----- |
| Supplementary Figure 1 | Positive control: BMI → T2D | `manuscript_figures.qmd` (`positive-control`) | `figures/ms-positive-control-1.pdf` | `results/pipeline/results.csv`, `results/pipeline/instruments_classical.csv` |
| Supplementary Figure 2 | Instrument strength (R², mean F) | `manuscript_figures.qmd` (`instrument-strength`) | `figures/ms-instrument-strength-1.pdf` | `results/pipeline/results.csv` |
| Supplementary Figure 3 | Concordance across 11 estimators | `manuscript_figures.qmd` (`method-concordance`) | `figures/ms-method-concordance-1.pdf` | `results/pipeline/results.csv` |
| Supplementary Figure 4 | SlopeHunter mechanism (A–D) | `slopehunter_figure.qmd` | `figures/supp-slopehunter-mechanism.pdf` | `results/slopehunter_figure/{fit_snps,b_refits,bmi_positions,bmi_genes}.csv` (computed once from the targets store; `-P refit:true` recomputes); `results/pipeline/{results,instruments_classical}.csv`; the Bellenguez b_SH fit |
| Supplementary Figure 5 | SlopeHunter sensitivity to b | `slopehunter_figure.qmd` | `figures/supp-slopehunter-sensitivity.pdf` | `results/slopehunter_figure/sensitivity_grid.csv`; `results/mvmr_lifespan/{life_slope,results}.csv`; `results/pipeline/bsh_null_checks.csv` |
| Supplementary Figure 6 | Survival collider by MVMR with lifespan (Zhang et al. 2026, step 2) | `mvmr_lifespan_indexevent.qmd` (`forest`) | `mvmr_lifespan_indexevent_files/figure-html/forest-1.png` — **PNG only; add `fig.path` and a PDF device before submission** | `results/mvmr_lifespan/*.csv` (`R/mvmr_lifespan.R`; `-P refit:true` refits from the targets store) |
| Supplementary Figure 7 | Fasting vs 2 h glucose: separability and the added-variable plot, both instrument sets | `manuscript_figures.qmd` (`fg-2h-validation`) | `figures/ms-fg-2h-validation-1.pdf` | computed in the chunk from `contrast/instrument_outcome_data.csv`, `contrast/step0_strength.csv` |
| Supplementary Figure 8 | Cluster estimates, primary vs replication outcome | `manuscript_figures.qmd` (`clusters-replication`) | `figures/ms-clusters-replication-1.pdf` | `contrast/cluster_results.csv` |
| Supplementary Figure 9 | Drug-target classes, primary vs replication outcome | `manuscript_figures.qmd` (`drug-targets-replication`) | `figures/ms-drug-targets-replication-1.pdf` | `contrast/drug_target_results.csv` |
| Supplementary Figure 10 | Combined Uffelmann 2026 (all cohorts) sensitivity analysis | `manuscript_figures.qmd` (`uffelmann-combined`) | `figures/ms-uffelmann-combined-1.pdf` | `contrast/{results,cluster_results,drug_target_results,instrument_outcome_data}.csv` |
| *(optional)* | KCNJ11/ABCC8 fine-mapping and colocalisation | `diabetes_ad_fine_mapping.qmd` | `figures/finemap-*-1.pdf` | **Reads the local `_targets` store on every render.** Save its inputs as committed CSVs before it goes in the paper |

Each `manuscript_figures.qmd` chunk also writes the data it plots to
`results/manuscript/source_data/<chunk label>.csv` (with a suffix for multi-panel chunks). Use these
for journals that ask for source data.

## Tables

| Table | Content | Script (chunk) | File |
| ----- | ------- | -------------- | ---- |
| Table 1 | Data sources (GWAS, accession, N, cohorts, covariate adjustment) | `manuscript_figures.qmd` (`table-data-sources`) | `results/manuscript/tables/data_sources.csv` (built from `pipeline/config.yml`) |
| Table 2 | Instrument strength: univariable, and conditional F for fasting vs 2 h glucose | `manuscript_figures.qmd` (`table-instrument-strength`) | `results/manuscript/tables/instrument_strength_univariable.csv`, `instrument_strength_fg_2h.csv` |
| Table 3 | Fasting vs 2 h glucose, primary outcome | `manuscript_figures.qmd` (`table-fg-2h`) | `results/manuscript/tables/fg_2h_primary_outcome.csv` |

| Supplementary Table | Content | Script (chunk) | File |
| ------------------- | ------- | -------------- | ---- |
| Supplementary Table 1 | All univariable MR estimates and diagnostics | targets pipeline (`csv_results`) | `results/pipeline/results.csv` (subset and relabel for the paper) |
| Supplementary Table 2 | Fasting vs 2 h glucose, every outcome and both scales | `manuscript_figures.qmd` (`table-fg-2h`) | `results/manuscript/tables/fg_2h_all_outcomes.csv` |
| Supplementary Table 3 | Cluster estimates, every outcome | `manuscript_figures.qmd` (`table-clusters`) | `results/manuscript/tables/cluster_results.csv`, `cluster_primary_outcome.csv` |
| Supplementary Table 4 | Drug-target classes, every outcome | `manuscript_figures.qmd` (`table-drug-targets`) | `results/manuscript/tables/drug_target_results.csv` |
| Supplementary Table 5 | KATP specificity steps, every outcome | `manuscript_figures.qmd` (`table-katp-specificity`) | `results/manuscript/tables/katp_specificity.csv` |
| Supplementary Table 6 | Instruments and their variant-level associations | targets pipeline | `results/pipeline/instruments_classical.csv`, `contrast/step0_instruments.csv`, `contrast/instrument_outcome_data.csv` |
| Supplementary Table 7 | SlopeHunter estimates at fixed b | targets pipeline (`csv_bsh_grid`) | `results/pipeline/bsh_grid.csv` |
| Supplementary File 1 | STROBE-MR checklist | *no script* | `Manuscripts/T2D_AD_MR/STROBE-MR_checklist.md` (working copy) → `STROBE-MR-checklist.docx` at submission |

## Where the data come from

```
GWAS summary statistics           pipeline/config.yml → datasets:  (Table 1; files in data/cache/, gitignored)
  │
  ├─ targets pipeline              _targets.R, pipeline/R/  (run on Great Lakes; locally only by target name)
  │    └─ results/pipeline/*.csv   history: results/pipeline/RUN_INFO.md · per-file stamps: results/pipeline/stamps/
  │
  ├─ Bellenguez b_SH fit           see "Survival-collider correction" below
  │
  └─ notebooks                     exposure_survival_screen.qmd → ad_longevity_arm.csv
                                   slopehunter_figure.qmd       → results/slopehunter_figure/
                                   mvmr_lifespan_indexevent.qmd → results/mvmr_lifespan/
        │
        └─ manuscript_figures.qmd  → figures/ms-*, results/manuscript/{tables,source_data}/
                                     results/manuscript/provenance.csv (each input, the last commit
                                     that changed it, uncommitted changes, pipeline stamp)
                                     results/manuscript/session_info.txt
```

**Survival-collider correction.** Every SlopeHunter-corrected estimate subtracts b_SH × (variant
effect on lifespan). The Bellenguez fit is `R/fits/bellenguez_lifespan_slopehunter_fits.csv` (set
`excl_APOE (PRIMARY)`: b_SH = −0.944, SE 0.070, 2,250 fitted SNPs), written by
`scripts/refit_bellenguez_bsh.R` from the GWAS Catalog harmonised Bellenguez file and the parental
lifespan GWAS. It is identical to full precision to the smoking-pipeline fit
(`smoking/results/indexevent_slopehunter_fits.csv`) that the results up to `7d9cad5` were computed
with. The pipeline reads it via `pipeline/config.yml` → `slopehunter: fits: ad_bellenguez`.

**Inputs from outside `Human Genetics/alzheimers/`:**

- `../calcium-cholesterol/raw_data/reference/w_hm3.snplist` (HapMap3 SNP list for MRAID and CAUSE)
- `../calcium-cholesterol/raw_data/reference/eur_w_ld_chr/` (LD scores for MR-APSS)
- `../calcium-cholesterol/R/robust_mr_helpers.R` (sourced by `_targets.R`)

**Outcome-blind design.** The fasting vs 2 h glucose instruments, cluster assignments and drug-target
classes were frozen before any AD data were read: `assignments/` with `MANIFEST.csv`, Step 0 results,
commit `38987ad` (2026-10-12). Outcome steps came after (`cccee6b`, 2026-10-13). Post-freeze
additions are labelled in the figures and in `STROBE-MR_checklist.md` (item 9b).

## Regenerating

From `Human Genetics/alzheimers/`:

```bash
quarto render manuscript_figures.qmd
```

For submission, refuse to render if any input has uncommitted changes:

```bash
quarto render manuscript_figures.qmd -P strict:true
```

If `quarto` is not on the PATH, use `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`.
Values still marked "TO CONFIRM" at their source are listed in the notebook's **Open items** section.
