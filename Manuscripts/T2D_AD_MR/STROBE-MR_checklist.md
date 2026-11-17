# STROBE-MR checklist — working copy

Tracks each STROBE-MR item (Skrivankova et al., *JAMA* 2021;326:1614–1621; *BMJ* 2021;375:n2233)
while the manuscript is written: where it is (or will be) addressed and which file in the repo backs
it. **Update the Section and Status columns as each part of `Manuscript.docx` is drafted.**

Item text here is a short paraphrase. For submission, fill in the official checklist from
<https://www.strobe-mr.org/>, adding page or line numbers, and save it here as
`STROBE-MR-checklist.docx` (README: Supplementary File 1).

Repo paths are relative to `Human Genetics/alzheimers/` unless they start with `Manuscripts/`.

**Status key:** `todo` = not yet in the manuscript · `partial` = some of it drafted · `done` = drafted
and checked against the repo · `n/a` = with the reason.

## Title and abstract

| Item | Asks for | Section | Status | Evidence in repo / notes |
|---|---|---|---|---|
| 1 | MR named as the design in the title and/or abstract | Title, Abstract | todo | — |

## Introduction

| Item | Asks for | Section | Status | Evidence in repo / notes |
|---|---|---|---|---|
| 2 | Background and rationale: the exposure, why a causal effect is plausible, why MR suits the question | Introduction | todo | — |
| 3 | Objectives, including pre-specified causal hypotheses; MR estimates causal effects only under assumptions | Introduction | todo | The fasting vs 2 h glucose decision rule was fixed before any AD data were read: `pipeline/config.yml` → `contrast: cond_f_min`, and `pipeline/README.md` ("Fasting vs 2 h glucose contrast"). Say which hypotheses were pre-specified and which came later (see 9b). |

## Methods

| Item | Asks for | Section | Status | Evidence in repo / notes |
|---|---|---|---|---|
| 4 | Key elements of the design early, ideally a table of data sources | Methods | todo | `results/manuscript/tables/data_sources.csv` (built from `pipeline/config.yml` by `manuscript_figures.qmd`) |
| 4a | Setting, locations and dates of each data source | Methods | todo | Per source GWAS paper; accessions in `data_sources.csv` |
| 4b | Eligibility and participant selection in each source | Methods | todo | Per source GWAS paper. Bellenguez includes family-history proxy cases; Kunkle and Uffelmann clinical-only are clinical |
| 4c | Measurement, quality control and selection of genetic variants | Methods: Instrument Selection | todo | `pipeline/README.md` ("Instrument sets", dataset QC); `config.yml` → `instruments:`, `contrast: instruments`, `contrast: drug_targets: cis`; APOE ± 500 kb exclusion (`contrast: apoe`) |
| 4d | How each exposure and outcome was assessed; diagnostic criteria for diseases | Methods | todo | Per source paper. Note the BMI adjustment of both Chen 2021 glucose GWAS (`config.yml` → `adjusted_for`) |
| 4e | Ethics approval and consent | Methods | todo | Public summary statistics only; cite each source's approvals. Suzuki 2024 (DIAGRAM click-through terms) and Uffelmann 2026 (CC BY-NC-SA 4.0) have use terms (`config.yml` comments) |
| 5 | The three core IV assumptions (relevance, independence, exclusion restriction), plus assumptions of any extra analyses | Methods | todo | DAGs: `mr_dags.qmd` (`iv-generic`, `collider`). SlopeHunter's own assumptions: `slopehunter_figure.qmd` |
| 6a | Scale, units and model for quantitative variables | Methods: Selection of MR Estimators | todo | Per mmol/L (MAGIC betas) vs per SD. **Per-SD estimates depend on the unconfirmed `sd_mmol`** (`manuscript_figures.qmd` → Open items). MR-APSS, CAUSE and MRAID are SD per SD; classical methods are per GWAS unit (`pipeline/README.md`, "Scales differ") |
| 6b | How variants were handled and weighted | Methods: Instrument Selection | todo | Thresholds, clumping and harmonisation: `pipeline/README.md`; Kunkle palindromes dropped (no allele frequencies) |
| 6c | The MR estimator(s), covariates, and whether the two samples were adjusted for the same covariates | Methods: Selection of MR Estimators | todo | Main MR-APSS, sensitivity CAUSE, MRBEE, IVW-MRE (residual SE floored at 1), etc.: `pipeline/README.md` ("Design"). MVMR and contrast regression: `pipeline/R/contrast_outcome.R`. **Covariate mismatch:** Chen glucose GWAS are BMI-adjusted; T2D and AD are not |
| 6d | Handling of missing data | Methods | todo | No LD-proxy substitution; missing instruments are dropped and logged: `results/pipeline/instrument_attrition.csv.gz`, `contrast/step0_attrition.csv` |
| 6e | Multiple testing | Methods | todo | **Decide:** which family of tests, and whether FDR is used (`SUMMARY.md` says the corrected arm "survives FDR"; the computation is not yet in the manuscript notebook) |
| 7 | How the assumptions were assessed or justified | Methods | todo | F and conditional F (`contrast/step0_strength.csv`), Steiger, Egger intercept, MR-PRESSO, Q/I² (`results/pipeline/results.csv`); BMI → T2D positive control; survival-collider screen (`exposure_survival_screen.qmd`) |
| 8 | Sensitivity and additional analyses (other estimators, replication, bias analyses, instrument validation) | Methods: Index Event Bias Adjustment | todo | b_SH uncertainty propagation and fixed-b grid (`results/pipeline/bsh_grid.csv`); MVMR with lifespan (`mvmr_lifespan_indexevent.qmd`); Kunkle and Uffelmann outcomes; p < 5×10⁻⁸ MR-APSS threshold; combined Uffelmann |
| 9a | Software, packages, versions and settings | Methods | todo | `results/manuscript/session_info.txt`; per-file package versions in `results/pipeline/stamps/` (pipeline runs from 2026-11-17 on) |
| 9b | Whether the protocol was pre-registered, and when and where | Methods | todo | Not externally registered. **Outcome-blind freeze in git:** `assignments/MANIFEST.csv` and Step 0 committed in `38987ad` (2026-10-12), before any outcome step (`cccee6b`, 2026-10-13). **Declare as post-freeze:** Uffelmann 2026 outcomes (2026-10-29), bromocriptine class (2026-10-15), the Xue → Suzuki cluster LD map, and the KATP specificity steps (post hoc) |

## Results

| Item | Asks for | Section | Status | Evidence in repo / notes |
|---|---|---|---|---|
| 10a | Numbers at each stage and reasons for exclusion; consider a flow diagram | Results | todo | GWAS N: `data_sources.csv`. Instrument counts and attrition: `contrast/step0_attrition.csv`, `instrument_attrition.csv.gz`, `results.csv` (`n_snp`) |
| 10b | Summary statistics for exposures and outcomes | Results | todo | From source papers. Glucose SDs for any per-SD statement (see 6a) |
| 10c | Heterogeneity across studies, where sources are meta-analyses | Results | todo | From source papers |
| 10d | Two-sample MR: similarity of exposure and outcome samples; overlap between them | Results | todo | All European ancestry. Overlap: `contrast/step0_overlap.csv`, `results/pipeline/associations.csv` (`overlap`), MR-APSS C12 (`results.csv` → `overlap_C12`). **MVP is in both Suzuki T2D and the Uffelmann new cohorts** |
| 11a | Variant–exposure and variant–outcome associations, on an interpretable scale | Results, Supplement | todo | `results/pipeline/instruments_classical.csv`, `contrast/step0_instruments.csv`, `contrast/instrument_outcome_data.csv` → supplementary table(s) |
| 11b | MR estimates with uncertainty, on an interpretable scale | Results | todo | `results/manuscript/tables/*.csv` (README table) |
| 11c | Absolute risk translation, if relevant | Results | todo | Probably n/a; decide |
| 11d | Plots, e.g. forest plots and variant scatterplots | Figures | todo | `figures/ms-*.pdf`, README figure table |
| 12a | Assessment of the validity of the assumptions | Results | todo | As in 7 |
| 12b | Additional statistics, e.g. heterogeneity (I², Q) | Results, Supplement | todo | `results.csv`: `Q`, `Q_p`, `I2`, `I2_GX`, `egger_intercept_p`, `presso_global_p`; cluster heterogeneity in `results/manuscript/tables/cluster_results.csv` |
| 13a | Sensitivity to violations of the assumptions | Results | todo | Method concordance (`figures/ms-method-concordance`); b_SH grid; MVMR with lifespan |
| 13b | Other sensitivity or additional analyses | Results | todo | Replication outcomes; combined Uffelmann (`figures/ms-uffelmann-combined`) |
| 13c | Direction of the causal relationship (e.g. bidirectional MR) | Results | todo | Reverse directions are in `results.csv` (`analysis: reverse_directions`); Steiger columns |
| 13d | Comparison with non-MR estimates, where relevant | Discussion | todo | Observational T2D/glucose–dementia literature |
| 13e | Further plots, e.g. leave-one-out | Supplement | todo | `results/pipeline/figures/<assoc_id>/leave_one_out.png` |

## Discussion

| Item | Asks for | Section | Status | Evidence in repo / notes |
|---|---|---|---|---|
| 14 | Key results, relative to the objectives | Discussion | todo | **The pre-specified primary (contrast regression, Suzuki set, MRBEE) is null:** OR 1.06 (0.95–1.20) per mmol/L (`results/manuscript/tables/fg_2h_primary_outcome.csv`). The Xue-set results are a sensitivity analysis |
| 15 | Limitations: IV assumptions, other biases, imprecision, with direction and size of bias | Discussion | todo | BMI-adjusted 2 h glucose GWAS (collider, `SUMMARY.md` §5); SlopeHunter Type I error (Zhang/Burgess 2026); conditional F < 10; MVP overlap; unconfirmed `sd_mmol` |
| 16a | Cautious overall interpretation in context of other studies | Discussion | todo | — |
| 16b | Mechanism; gene–environment equivalence; careful causal language | Discussion | todo | KATP and beta-cell findings; fine-mapping (`diabetes_ad_fine_mapping.qmd`) |
| 16c | Clinical or policy relevance; what it implies for intervention effect sizes | Discussion | todo | Drug-target results (sulfonylureas, glucokinase) |
| 17 | Generalisability across populations, timing and levels of exposure | Discussion | todo | European ancestry only; lifelong genetic exposure |

## Other information

| Item | Asks for | Section | Status | Evidence in repo / notes |
|---|---|---|---|---|
| 18 | Funding for this study and for the source data | Acknowledgements and Funding | todo | — |
| 19 | Data and code: where they can be accessed; code to reproduce the results | Data and code availability (add) | todo | Repo plus a tagged release at submission; `Manuscripts/T2D_AD_MR/README.md` maps every figure and table to its code. GWAS accessions in `data_sources.csv`. Suzuki and Uffelmann files cannot be reposted |
| 20 | Conflicts of interest | Declarations | todo | — |
