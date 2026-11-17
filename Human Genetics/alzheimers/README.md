# Alzheimer's survival-collider analysis

Mendelian-randomization workup of the **survival collider** in Alzheimer's-disease
(AD) MR. Because an AD outcome GWAS is measured only in people who *survived to be
assessed*, the conditioning variable is **survival `S`**, not AD incidence.
Conditioning on `S` (studying survivors) opens the path

```
   G_exposure ──▶ S ◀── AD        (AD → death is the known edge)
                  ▲
                  └── competing risks (IHD, stroke, cancer) sharing AD etiology
```

and induces a spurious `G_exposure`–AD association (index-event / survival bias).
This folder screens each exposure for that bias, then corrects the exposure→AD
estimates for it.

## Pipeline

| Step | Notebook | Question |
|------|----------|----------|
| **1. Screen** | [`exposure_survival_screen.qmd`](exposure_survival_screen.qmd) | Does each exposure have an edge into the survival collider (`E → S`)? Establishes the `AD → S` arm, screens `E → S`, and nominates the induced-bias direction with a Bayesian ROPE. |
| **2. Correct** | [`exposure_ad_indexevent.qmd`](exposure_ad_indexevent.qmd) | Naive `E → AD` (run through the per-association SOP-validated model set in the config — IVW-MRE, MR-PRESSO, MR-RAPS, …), then SlopeHunter-corrected `E → AD` = naive − `b_SH·(E→lifespan)`. Flags each exposure ABOLISHED / ROBUST / UNMASKED / REVERSED. |

The correction reuses the exact approach validated in the smoking→AD analysis
(`smoking/04e_overlap_indexevent.qmd`). `b_SH` is the SlopeHunter selection slope
relating SNP effects on AD to SNP effects on lifespan; it is a property of the
**(AD GWAS × lifespan GWAS)** pair and is exposure-independent, so it is fit once
per AD GWAS and applied to every exposure. Its uncertainty is carried into every
corrected estimate: Step 2 adds `(E→lifespan)²·se(b_SH)²`, and the pipeline's
SlopeHunter arm propagates it per estimate (see `pipeline/README.md`, "Survival collider").

## Files

```
exposure_survival_screen.qmd        Step 1 — screen (self-contained; OpenGWAS API)
exposure_ad_indexevent.qmd          Step 2 — naive + SlopeHunter correction
R/fit_selection_slope.R             SlopeHunter prep + fit + linear-correction engine
mvmr_lifespan_indexevent.qmd        Cross-check — MVMR with lifespan as a second exposure (Zhang et al. 2026, step 2)
R/mvmr_lifespan.R                   its engine; writes results/mvmr_lifespan/
slopehunter_figure.qmd              Supplementary figures: SlopeHunter mechanism (A–D) and sensitivity to b → figures/supp-slopehunter-*
manuscript_figures.qmd              Greyscale manuscript figures and tables → figures/ms-*, results/manuscript/ (key: Manuscripts/T2D_AD_MR/README.md)
t2d_ad_figures.qmd                  Colour talk versions of the same panels → figures/talk-*
scripts/refit_bellenguez_bsh.R      Bellenguez b_SH fit → R/fits/bellenguez_lifespan_slopehunter_fits.csv (identical to the smoking-pipeline fit)
scripts/fetch_kunkle_sumstats.sh    Turnkey fetch/normalize of Kunkle b_SH inputs
docs/MR_visual_inspection.mediawiki SOP section (scatter/funnel + decision matrix) for the lab wiki
```

Outputs written on render: `exposure_survival_screen.csv`,
`exposure_survival_verdict.csv`, `ad_longevity_arm.csv`,
`exposure_ad_naive.csv` (cached naive estimates), and
`exposure_ad_indexevent_corrected.csv`.

## Requirements

- **R** with `TwoSampleMR`, `ieugwasr`, `dplyr`, `tidyr`, `purrr`, `readr`,
  `ggplot2`, `ggdag`; plus `SlopeHunter` and `data.table` for the correction.
  Optionally `MRPRESSO` and `mr.raps` — only needed if an exposure's `naive_model`
  is `MR-PRESSO` or `MR-RAPS` (absent → the model falls back to IVW-MRE).
- **An OpenGWAS JWT** for the API steps. Get one at <https://api.opengwas.io>, add
  `OPENGWAS_JWT=<token>` to `~/.Renviron`, and restart R. Treat it like a password.
- **`plink2`** on your PATH — only for the LD-clumping step when fitting a new
  `b_SH` (Kunkle below). `quarto` to render.

## How to run

Run everything from **this folder** (script paths are relative to it):

```bash
cd "Human Genetics/alzheimers"
```

### Step 1 — screen

```bash
quarto render exposure_survival_screen.qmd
```

### Step 2 — naive + corrected exposure → AD

The naive `E → AD` step needs a valid OpenGWAS token (it caches to
`exposure_ad_naive.csv` on first success). The Bellenguez `b_SH` is read from
`R/fits/bellenguez_lifespan_slopehunter_fits.csv` (`Rscript scripts/refit_bellenguez_bsh.R`;
identical to the fit first made in the smoking pipeline), so rendering needs no local
sumstats. Kunkle rows show "run the fit for this AD GWAS" until you do the fit below:

```bash
quarto render exposure_ad_indexevent.qmd
```

### Fitting `b_SH` for Kunkle (or any new AD GWAS)

`b_SH` needs full genome-wide summary statistics (the OpenGWAS API only serves
clumped tophits), so it is fit outside the notebook. For Kunkle 2019 IGAP the
fetch script is turnkey — it downloads the EBI-hosted IGAP Stage-1 file, normalizes
its columns to the harmonised schema, and reuses your smoking cache for the shared
lifespan sumstats + LD panel:

```bash
bash scripts/fetch_kunkle_sumstats.sh
```

The script prints the remaining commands; they are (run from this folder):

```bash
# 1. prep + build the plink clump input
Rscript -e 'source("R/fit_selection_slope.R"); m <- prep_selection_merge("data/cache/sumstats/AD_kunkle_harmonised.tsv.gz","data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz"); data.table::fwrite(m[,.(ID=SNP,P=p_life)],"data/cache/sumstats/kunkle_clump_input.tsv",sep="\t")'
```

```bash
# 2. LD-clump to an independent SNP set
cd data/cache && plink2 --bfile EUR --clump sumstats/kunkle_clump_input.tsv --clump-p1 1e-3 --clump-p2 0.01 --clump-r2 0.01 --clump-kb 1000 --out sumstats/kunkle_clumped && cd ../..
```

```bash
# 3. fit SlopeHunter -> writes R/fits/kunkle_lifespan_slopehunter_fits.csv
Rscript -e 'source("R/fit_selection_slope.R"); dir.create("R/fits",showWarnings=FALSE,recursive=TRUE); print(build_selection_slope(ad_file="data/cache/sumstats/AD_kunkle_harmonised.tsv.gz", life_file="data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz", clumped_snps="data/cache/sumstats/kunkle_clumped.clumps", out_csv="R/fits/kunkle_lifespan_slopehunter_fits.csv"))'
```

Then re-render `exposure_ad_indexevent.qmd` and the Kunkle rows populate.

## Reading the correction

- **ABOLISHED** — significant naive, but the corrected CI covers 0; the effect is
  explained by the survival collider. Do not treat it as a causal/druggable signal.
- **ROBUST** — significant before *and* after correction; the collider does not
  account for it.
- **UNMASKED** — null naive but significant after correction; the collider was
  *suppressing* a real effect. A null naive never means "no effect."
- **REVERSED** — correction flips the sign; the induced bias dominated the naive
  estimate.

The Step-1 `induced_bias` direction predicts *which way* the naive estimate is
pulled; Step 2 quantifies *how much*. They should agree in sign. **UNMASKED and
REVERSED results are only as reliable as `b_SH` and `E→lifespan`** — a mis-oriented
selection slope can manufacture a "revealed" effect, so verify those inputs (allele
orientation, `b_SH` sign) before believing one.

**Caveats.** `b_SH` is per AD GWAS (fit each; don't reuse one across ascertainments).
The correction cannot separate collider-induced association from a genuine effect
collinear with the mortality axis — strong evidence, not proof. Kin-proxy parental
lifespan halves per-allele effects, so the subtracted term is conservative.

## Known issues

- **Kunkle's `b_SH` is not identified, so it is never applied.** The fit
  (`R/fits/kunkle_lifespan_slopehunter_fits.csv`, committed) returns `+1.10`,
  opposite in sign to Bellenguez's `−0.944`. `exposure_ad_indexevent.qmd` ("Why the
  Kunkle `b_SH` is not used") ruled out allele inversion and a sign-convention
  mismatch. On the fitted SNPs there is no raw lifespan–AD relationship (r = −0.0008,
  sign agreement 0.499), so SlopeHunter is partitioning noise. The Uffelmann 2026
  new-cohorts and clinical-only files are not identified either (raw r = 0.034, 0.047).
  The combined Uffelmann file is weak (r = 0.115), and its `b_SH` sits inside the
  permutation null. The Kunkle and Uffelmann columns are reported uncorrected;
  `pipeline/config.yml` marks all four fits `trusted: false`.
  Run `bsh_identifiability()` on any new AD GWAS before trusting its `b_SH`.
