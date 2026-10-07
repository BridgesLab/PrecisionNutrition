# Univariable robust-MR pipeline (`targets`)

Re-estimates every exposure → mediator, exposure → AD and exposure → lifespan association with a
pre-specified **primary** method that allows correlated pleiotropy, plus the full classical suite
and diagnostics, from full summary statistics, with every SNP loss logged. It also renders the
existing Step 1, Step 2 and mediation notebooks in order.

Definition: [`../_targets.R`](../_targets.R). Settings: [`config.yml`](config.yml).

## Design

| Piece | Rule |
|---|---|
| **Main method** | **MR-APSS** for every association. It models correlated pleiotropy (background genetic correlation Ω) *and* sample overlap (LDSC intercepts C), so one method covers UKB-overlapping and non-overlapping pairs alike. Instruments p < 5×10⁻⁵ with the winner's-curse correction, clumped r² < 0.001 / 1 Mb. |
| **Sensitivity** | **CAUSE**, also category 1, for every association. Its 95% credible interval for γ and P(γ < 0) are the decision rule; the ELPD test is reported alongside. |
| **Category-2 comparison** | **MRBEE** and **MR-RAPS** (assume InSIDE; more powerful), always reported. They are flagged unstable only on process grounds (fit error, no finite SE, or RAPS's over-dispersion optimiser failing), never on their estimates. See `R/stability.R`. |
| **Threshold sensitivity** | **MR-APSS at p < 5×10⁻⁸** (config `mrapss_sensitivity`), separating "the robust model discounts the signal" from "the 5×10⁻⁵ instrument set dilutes it". |
| **Latent confounding** | CAUSE **q·η** (joint posterior draws) and q·η/(γ + q·η); MR-APSS **background slope Ω_xy/Ω_xx** with a 200-block jackknife interval, plus the sample-structure slope C_xy/C_xx. See `R/latent.R`, which logs the exact package object paths used. |
| **Reverse directions** | Every non-classical-only pair is also run reversed (raw arm), because reverse causation can load onto the latent terms. |
| **MRAID** (extra sensitivity) | Run only where (1) the two GWAS share no cohort (MRAID assumes independent samples) and (2) there are at most 1,000 candidate instruments (`analysis: mraid_rule`). Pairs over the cap appear as "not run" with the count. Reason: for lean mass and fat-free mass (3.5–5.3k correlated candidates, 1000G LD) MRAID's estimates were **unreliable**: effects > 1 SD/SD, and opposite signs for two fasting-glucose GWAS. With 25–300 candidates (glucose traits, T2D) it behaved well. |
| **Classical suite** | IVW-MRE, IVW-FE, MR-Egger (+ intercept), weighted median, weighted mode, MR-RAPS, MR-PRESSO, MRBEE, on the classical instrument set. |
| **Diagnostics** (every classical row) | Q, Rücker's Q′, I², I²_GX, R² (z-based), mean F, total F, Steiger (z-based R², approximate for binary traits). |
| **Survival collider** | Outcomes whose Step 2 b_SH is trusted (Bellenguez) get a second, SlopeHunter-adjusted arm, built genome-wide as β − b_SH·β_lifespan. Kunkle is never adjusted. The Bellenguez b_SH is checked against its permutation null (`bsh_null_checks.csv`). |
| **Step 2 gap** | `step2_lifespan_gap_fill.csv` gives E → lifespan for every exposure and Step 2's corrected E → AD, which was missing for lean mass. |

### Instrument sets

| Set | Used by | Selection |
|---|---|---|
| classical | classical suite | p < 5×10⁻⁸; clump r² < 0.001, 10 Mb (prior univariable criteria) |
| mraid | MRAID | HapMap3; p < 5×10⁻⁸; clump r² < 0.5, 1 Mb (stepping down to 0.3 / 0.2 / 0.1 only if still above MRAID's 10,000-SNP cap) |
| cause | CAUSE | HapMap3 genome-wide for nuisance parameters; p < 10⁻³, prune r² < 0.01, 10 Mb |

Before any of these, every dataset goes through the same QC: rsID in the 1000G EUR panel,
biallelic SNV, valid beta/SE/p, MAF ≥ 0.01, one row per rsID.

### Deviations from the prior univariable pipeline (all forced or deliberate)

- **Clumping runs locally** (plink2, same 1000G EUR panel) instead of on the OpenGWAS LD server,
  which has returned 502 since 2026-10-01. Checked: lean mass → Kunkle by IVW-MRE gives
  −0.126 (SE 0.038, 559 SNPs) against Step 2's −0.130 (SE 0.039, 528 SNPs).
- **No LD-proxy substitution.** Instruments missing from an outcome are dropped and logged.
  Proxy lookup also depended on the down LD service.
- **IVW-MRE floors the residual SE at 1**, the standard multiplicative random-effects model
  (matches `MendelianRandomization::mr_ivw(model = "random")`). TwoSampleMR 0.7.5's `mr_ivw_mre()`,
  which Step 2 uses, does not floor: with under-dispersed instruments (Q < df) its SE is smaller
  than fixed effects (about half, in a test with dispersion 0.6). MR-PRESSO's SE is unfloored by
  that method's own design and is left as is.
- **p-value underflow:** p-values below the double range (Kunkle APOE, 10⁻⁸⁸¹) are floored rather
  than dropped, and counted in the dataset QC log.

### Data notes

- **Bellenguez (`ebi-a-GCST90027158`) does not contain rs429358 or rs7412**, the two SNPs that
  define APOE ε2/ε3/ε4. The rest of the APOE region is present. This doesn't affect SlopeHunter
  (APOE is excluded from the fit by design); any exposure instrument at those SNPs is logged as
  missing in the outcome.
- **Kunkle has no allele frequencies**, so palindromic instruments are dropped as ambiguous under
  the harmonisation rule (78 of 637 for lean mass), exactly as in the prior API-based analyses.

## Outputs (`results/pipeline/`)

| File | Contents |
|---|---|
| `results.csv` | one row per association × method; `is_primary` marks the pre-specified estimate; diagnostics on classical rows | Columns for filtering overlap-naive methods: `overlap` (shared cohorts from the config), `overlap_C12` (MR-APSS's empirical cross-trait LDSC intercept) and `method_handles_overlap` (TRUE for MR-APSS, CAUSE, MRBEE).
| `instruments_primary.csv` | the SNPs each association's primary method used (alleles, effects on both traits) |
| `instrument_attrition.csv.gz` | every step from file → final set, with counts and the SNPs removed at each step |
| `associations.csv` | association list, cohort overlap, primary method |
| `datasets.csv` | datasets, sample sizes, cohorts, covariate adjustments (collider-risk flags) |
| `step2_lifespan_gap_fill.csv`, `bsh_null_checks.csv` | survival-collider outputs above |
| `instruments_classical.csv` | every association's classical instrument set (for custom scatter/LOO/funnel plots) |
| `mrapss_background.csv` | MR-APSS background per association: C12 (sample overlap), rg (genetic correlation) |
| `latent_confounding.csv` | one row per pair × method (CAUSE, MR-APSS): causal estimate, latent term with interval, fraction shared, CAUSE q / ΔELPD, MR-APSS sample-structure slope, flags |
| `cause_posteriors.rds` | CAUSE grid posteriors (causal and sharing) per pair, so latent terms can be recomputed without re-running CAUSE |
| `figures/<assoc_id>/` | `scatter_classical.png`, `funnel.png` (Wald ratios, IVW line, pseudo-95% cone), `leave_one_out.png` + `.csv`, `scatter_mrapss.png` |

**Report:** [`../robust_mr_report.qmd`](../robust_mr_report.qmd) reads only these files: datasets, a
main/sensitivity/category-2 table, forest plots split by scale, overlap, survival collider, and
tabbed diagnostics per association. Pipeline target `nb_report` (local only); it can also be
rendered on its own.

**Scales differ.** Classical estimates are per unit of each GWAS's beta (log-OR for binary
traits). MRAID and CAUSE are per SD of exposure on per SD of outcome. Compare across method
families on sign and significance, not magnitude.

## Running

```r
# local: everything except CAUSE (needs loo < 2.10) and the API-dependent notebooks
targets::tar_make(names = !starts_with("cause_") & !starts_with("nb_"))
targets::tar_read(results)
```

**Great Lakes** runs everything, including CAUSE. Copy this folder plus the pieces it reads from
outside it (`../calcium-cholesterol/R/robust_mr_helpers.R`, the HapMap3 list), then:

```bash
mkdir -p ~/logs logs && sbatch pipeline/run_greatlakes.sbatch
```

One-time package setup on the cluster, after `module load R gcc`:

```bash
Rscript -e 'install.packages(c("targets","tarchetypes","crew","crew.cluster","qs2","quarto","data.table","tidyverse","yaml","remotes","R.utils","mvtnorm"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages(c("MRPRESSO","mr.raps"), repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))'
Rscript -e 'remotes::install_github(c("MRCIEU/TwoSampleMR","stephenslab/mixsqp","stephens999/ashr","jean997/cause","yuanzhongshang/MRAID","noahlorinczcomi/MRBEE","Osmahmoud/SlopeHunter","YangLabHKUST/MR-APSS"))'
Rscript -e 'remotes::install_version("loo", version = "2.4.1", repos = "https://cloud.r-project.org")'
```

**plink2 on Great Lakes:** the cluster's plink2 modules (newest `plink/2.4a4`, i.e. 2.00a4 from
Jan 2023) have no `--clump`, which CAUSE's LD pruning needs. Install the current release in
`~/bin` once; the sbatch and the worker jobs use it via `PLINK2_BIN`:

```bash
mkdir -p ~/bin && cd ~/bin && wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20261001.zip && unzip -o plink2_linux_x86_64_20261001.zip plink2 && rm plink2_linux_x86_64_20261001.zip && ./plink2 --version
```

`mvtnorm` is in that list because MR-APSS calls `mvtnorm::dmvnorm()` internally without declaring
it as a dependency, so installing MR-APSS doesn't pull it in. Only some fits reach that code path.

Local note: MRAID was built against Homebrew's gfortran (`R_MAKEVARS_USER` pointing `FLIBS` at
`/opt/homebrew/opt/gcc/lib/gcc/current`) because CRAN's gfortran is not installed on this Mac.
