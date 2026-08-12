# ROBUST_MR.md — overlap- and pleiotropy-robust MR for cholesterol → BMD

Companion to `ANALYSIS.md` §3. Covers the arm implemented in
`robust_mr_prep.qmd`, `robust_mr_apss.qmd`, `robust_mr_cause.qmd`, and
`robust_mr_summary.qmd`.

---

## The problem this arm exists to solve

The heel BMD result (all-cholesterol → eBMD, β ≈ −0.051, p = 1.4e-4) is estimated
with a UK Biobank exposure against a UK Biobank outcome. That creates three
exposures to bias that IVW, Egger, weighted median and mode-based estimators
cannot address, because none of them model the joint distribution of the two GWAS:

1. **Sample overlap.** With near-complete overlap the analysis is effectively
   one-sample MR, and estimates are pulled toward the observational association.
2. **Correlated horizontal pleiotropy.** Lipids and bone share a substantial
   polygenic background (adiposity, growth, inflammation). A shared factor
   produces a proportional, low-heterogeneity pattern that Egger and weighted
   median read as a clean causal effect.
3. **Weak instrument bias / winner's curse.** Relaxing the p-value threshold to
   gain instruments normally imports selection bias.

## Method roster and why

| Method | Overlap | Correlated pleiotropy | Weak instruments | Status |
|---|---|---|---|---|
| **MR-APSS** | ✅ 2×2 `C` matrix from bivariate LDSC intercepts | ✅ background model `Ω` | ✅ 5e-5 threshold **with** explicit selection-bias correction | **Implemented** |
| **CAUSE** | ✅ nuisance parameter `ρ` estimated genome-wide | ✅ shared factor `U` (`η`, `q`) — its core purpose | ⚪ uses p<1e-3 variants, modelled not filtered | **Implemented** |
| **MRAID** | ❌ strict two-sample; assumes independent cohorts | ✅ spike-slab, correlated + uncorrelated | ✅ fine-mapping-style automated selection | **Deferred — see below** |

### Why MRAID is not here

MRAID (Yuan et al., *Science Advances* 2022) is a good method for the pleiotropy
and instrument-selection problems, but its likelihood is built on the standard
two-sample assumption that the exposure and outcome GWAS come from independent
samples. Running it on UK Biobank LDL against UK Biobank eBMD would violate the one
assumption it has no defence against, and the estimate would be biased toward the
observational association — the same failure the other two methods are here to
prevent.

**Planned extension.** MRAID becomes appropriate under a genuinely non-overlapping
design: **UK Biobank exposure → MGI / BioVU outcome**. That is the same
architecture already used for serum calcium in `mr-tc-calcium.qmd` and
`mr-ldlc-calcium.qmd`, so the instruments and harmonisation code largely exist. It
needs an MGI or BioVU BMD/fracture phenotype, which is the open question.

Recorded here as a deliberate exclusion with a stated remedy, not a gap.

## Study design — two arms, on purpose

| Arm | Exposure | N | UK Biobank? | Purpose |
|---|---|---|---|---|
| `overlapping` | LDL-C, GLGC 2021 (Graham, EUR) | ≈1.32M | Yes | Maximum power; the arm the correction is *for* |
| `overlap-free` | LDL-C, GLGC 2013 (Willer, `jointGwasMc_LDL.txt.gz`) | ≈173k | No — predates UK Biobank | Control; overlap is zero by construction |

Outcome for both: **heel eBMD, Morris 2019, UK Biobank, N = 426,824.**

Running both arms through the same pipeline turns "is overlap biasing us?" from an
argument into a measurement. Specifically:

- MR-APSS reports `C12`, the cross-trait LDSC intercept, which is a **direct
  empirical estimate of sample overlap** (× phenotypic correlation).
- CAUSE reports `ρ`, the null correlation of test statistics, which measures the
  same thing through a different model.
- MR-APSS is additionally fitted with `C = I` (correction off), so the shift
  between the two fits *is* the overlap bias in outcome units.

> **Open item.** `ANALYSIS.md` §3 states that the GLGC total-cholesterol exposure
> (`ebi-a-GCST90025953`) is non-overlapping with UK Biobank. UK Biobank does
> contribute to the Graham 2021 European meta-analysis, so this claim needs
> checking against Supplementary Table 1 of that paper. The `C12` estimate from
> this pipeline settles it empirically either way.

## Scope — what this arm does *not* touch

CAUSE and MR-APSS are genome-wide polygenic methods. They **cannot** be applied to
the cis-instrument designs in `drug_target_mr_bmd.qmd` (HMGCR, PCSK9, NPC1L1) —
those use a handful of SNPs in a single gene window and have a different
identifying assumption. This arm tests the *genome-wide circulating LDL-C → BMD*
claim only.

The interaction between the two is the interesting part. If genome-wide LDL → BMD
weakens under these methods while the HMGCR cis effect (β ≈ −0.115) holds, that
strengthens the paper's central argument that the mechanism is mevalonate-specific
rather than LDL-mediated.

## Running it

Everything is local except CAUSE. Only ~150 MB of harmonised cache ever needs to
move to the cluster, so the 12 GB of raw downloads stay on your machine.

### Step 0 — dependencies (local, once)

```bash
Rscript -e 'source("R/robust_mr_helpers.R"); check_robust_mr_deps(install = TRUE)'
```

CAUSE and MR-APSS are GitHub-only and pull `mixsqp`/`ashr`, which need a working
compiler. Confirm before going further:

```bash
Rscript -e 'library(cause); library(MRAPSS); cat("both load\n")'
```

### Step 1 — data (local, ~20–60 min, ~12 GB)

```bash
bash scripts/fetch_robust_mr_data.sh
```

Read what it prints at the end. Anything it could not fetch is listed with a manual
link — GEFOS and the Broad LD score host both move periodically.

### Step 2 — prep (local, ~20 min) — **two gates here**

```bash
quarto render robust_mr_prep.qmd
```

> **Gate 1.** The first chunk prints the real headers of every input file. Compare
> them against the `cols:` blocks in `config_robust_mr.yml` and fix any mismatch
> before continuing. A wrong effect-allele mapping is a silent sign flip.
> (The GLGC 2013 arm is already verified; GLGC 2021 and GEFOS are not.)
>
> **Gate 2.** The last chunk runs a plain IVW on the harmonised data. It must land
> near the −0.051 in `ANALYSIS.md` §3, with mean F well above 10. If it does not,
> the harmonisation is wrong and nothing downstream is interpretable. Stop and
> debug here, not later.

Also glance at the mean chi-square in the datasets table: below ~1.02 and the
LDSC-based background model will be unstable.

### Step 3 — MR-APSS (local, minutes)

```bash
quarto render robust_mr_apss.qmd
```

Do this before CAUSE. It is fast, and its `C12` estimate tells you immediately
whether the overlap problem is real — which is worth knowing before spending
compute on CAUSE.

### Step 4 — CAUSE (Great Lakes, ~1–2 h per arm, or overnight locally)

Fill in `#SBATCH --account=` first. Keep the project on `/scratch`, not home —
Great Lakes home directories have an 80 GB quota.

**The sbatch file alone is not enough.** Six things have to be there:

| Path | Why |
|---|---|
| `scripts/cause_greatlakes.sbatch` | the wrapper |
| `scripts/run_cause_greatlakes.R` | the actual job |
| `R/robust_mr_helpers.R` | sourced by the job for `load_robust_cfg()` and `clump_local()` |
| `config_robust_mr.yml` | paths, seed, pruning thresholds |
| `raw_data/robust_mr/cache/*_cause.rds` | the harmonised inputs from step 2 (~47 MB) |
| `alzheimers/data/cache/EUR.{bed,bim,fam}` | LD panel for pruning (~1.35 GB) |

None of the raw GWAS downloads are needed; `robust_mr_prep.qmd` already distilled
them. The LD panel is the bulk of the transfer.

> **Note on the LD panel path.** Locally the panel lives in a *sibling* directory
> (`../alzheimers/data/cache/EUR`), shared with the AD analyses, because this
> project moved into `calcium-cholesterol/`. On the cluster it was uploaded flat,
> *under* the project root. `resolve_bfile()` tries the configured path, then `../`
> and `../../`, so the single `plink_bfile` entry in `config_robust_mr.yml` works
> in both places. Absolute paths also pass through unchanged.

```bash
REMOTE=greatlakes.arc-ts.umich.edu:/nfs/turbo/sph-davebrid/GWAS_Calcium/robust_mr

rsync -avR \
  scripts/cause_greatlakes.sbatch \
  scripts/run_cause_greatlakes.R \
  R/robust_mr_helpers.R \
  config_robust_mr.yml \
  raw_data/robust_mr/cache/ \
  "${REMOTE}/"

# the LD panel is a sibling of this directory, so send it separately
rsync -av ../alzheimers/data/cache/EUR.{bed,bim,fam} \
  "${REMOTE}/alzheimers/data/cache/"
```

`-R` matters: it preserves the relative paths so `config_robust_mr.yml` finds
everything. Without it the files flatten into one directory.

One-time setup on the cluster:

```bash
ssh greatlakes.arc-ts.umich.edu
cd /nfs/turbo/sph-davebrid/GWAS_Calcium/robust_mr
mkdir -p ~/logs               # slurm opens --output BEFORE the script runs; without
                              # this the job dies with no error message

module load R gcc Bioinformatics plink/1.9
Rscript -e 'install.packages(c("remotes","tidyverse","data.table","yaml","ieugwasr"), repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_github(c("stephenslab/mixsqp","stephens999/ashr","jean997/cause"))'
# CAUSE 1.2.0 indexes loo_compare() positionally; loo >= 2.3 broke that.
Rscript -e 'remotes::install_version("loo", version="2.4.1", repos="https://cloud.r-project.org")'
```

`MRAPSS` is **not** needed on the cluster — MR-APSS runs locally in step 3.

```bash
sbatch --export=ALL,ARM=overlapping  scripts/cause_greatlakes.sbatch
sbatch --export=ALL,ARM=overlap-free scripts/cause_greatlakes.sbatch
squeue -u $USER
```

The job checks for its own inputs and exits immediately with a `MISSING:` line if
the transfer was incomplete, rather than failing 40 minutes in.

Both arms are independent, so submit them together. Then pull the fits back:

```bash
rsync -av "${REMOTE}/raw_data/robust_mr/cache/cause_fit_*.rds" \
  raw_data/robust_mr/cache/
quarto render robust_mr_cause.qmd    # detects the .rds, skips the slow chunks
```

To skip the cluster entirely, just render `robust_mr_cause.qmd` locally and leave
it running — it computes the fits itself if no cache is present.

### Step 5 — synthesis (local, seconds)

```bash
quarto render robust_mr_summary.qmd
```

Then fill the `[FILL: ...]` placeholders in the suggested `ANALYSIS.md` §3b text at
the bottom of that file.

### Compute profile

Sized for the HapMap3-restricted variant set (~1M SNPs) that `robust_mr_prep.qmd`
produces. The HM3 restriction is required by MR-APSS's LDSC step and happens to be
exactly the variant count CAUSE wants for `est_cause_params()`, so both methods run
on the same, comparatively small, input.

| Step | Time | Memory | Where |
|---|---|---|---|
| Download | 20–60 min | — | local (~12 GB disk) |
| QC + harmonisation (`prep`) | 10–20 min | 8–16 GB | local |
| MR-APSS `est_paras` (bivariate LDSC) | 2–5 min | 8 GB | local |
| MR-APSS fit ×4 | < 5 min | 4 GB | local |
| CAUSE `est_cause_params` + fit, per arm | 1–2 h | 16–32 GB | either |
| Synthesis | seconds | — | local |

Nothing here strictly requires a cluster. CAUSE is the only piece worth sending to
Great Lakes, and mainly so the two arms run in parallel instead of tying up your
laptop for an evening. `scripts/cause_greatlakes.sbatch` needs your slurm account
filled in (`#SBATCH --account=`).

## Files

| File | Role |
|---|---|
| `config_robust_mr.yml` | Paths, column mappings, thresholds. **Column mappings are unverified guesses — check them.** |
| `R/robust_mr_helpers.R` | Reading, QC, format conversion, clumping, result tidiers, forest plot |
| `scripts/fetch_robust_mr_data.sh` | Downloads with fallback URLs; reports what it could not get |
| `scripts/run_cause_greatlakes.R` | Headless CAUSE fit for one arm |
| `scripts/cause_greatlakes.sbatch` | Slurm wrapper |
| `robust_mr_prep.qmd` | Header inspection, QC, harmonisation, IVW sanity anchor |
| `robust_mr_apss.qmd` | MR-APSS, both arms, C estimated vs C = I |
| `robust_mr_cause.qmd` | CAUSE, both arms |
| `robust_mr_summary.qmd` | Reconciliation, forest plot, pre-registered interpretation grid |

## Known soft spots

- **Column mappings in `config_robust_mr.yml` are best guesses.** GLGC and GEFOS
  headers have changed between releases. `robust_mr_prep.qmd` prints the real
  headers first; correct the config before trusting anything. A wrong effect-allele
  column is a silent sign flip.
- **Download URLs** for GEFOS, the Broad LD scores, and GLGC have all moved in
  recent years. The fetch script tries fallbacks and reports failures rather than
  proceeding quietly.
- ~~`MRAPSS` result field names~~ — **verified 2026-08-10** against the installed
  package using the bundled BMI→T2D example. The object carries `beta`, `beta.se`,
  `pvalue` (**character**, not numeric — `tidy_apss()` coerces), `tau.sq`,
  `sigma.sq`, `pi0`, `IVsignal.sum`, `Threshold`, `MRdat`, `post$Pi`. There is no
  instrument-count field; it is `nrow(MRdat)`. `MRAPSS::clump()` exists and its
  signature matches the call in `robust_mr_apss.qmd`.
- The IVW sanity check in `robust_mr_prep.qmd` must reproduce something close to
  the −0.051 in `ANALYSIS.md`. If it does not, the harmonisation is wrong and
  nothing downstream means anything.

## References

- Morrison J, Knoblauch N, Marcus JH, Stephens M, He X. Mendelian randomization
  accounting for correlated and uncorrelated pleiotropy using genome-wide summary
  statistics. *Nat Genet* 2020;52:740–7. — CAUSE
- Hu X, Zhao J, Lin Z, Wang Y, Peng H, Zhao H, Wan X, Yang C. Mendelian
  randomization for causal inference accounting for pleiotropy and sample structure
  using genome-wide summary statistics. *PNAS* 2022;119(28):e2106858119. — MR-APSS
- Yuan Z, Liu L, Guo P, Yan R, Xue F, Zhou X. Likelihood-based Mendelian
  randomization analysis with automated instrument selection and horizontal
  pleiotropic modeling. *Sci Adv* 2022;8:eabl5744. — MRAID
- Graham SE, et al. The power of genetic diversity in genome-wide association
  studies of lipids. *Nature* 2021;600:675–9. — GLGC 2021
- Willer CJ, et al. Discovery and refinement of loci associated with lipid levels.
  *Nat Genet* 2013;45:1274–83. — GLGC 2013
- Morris JA, et al. An atlas of genetic influences on osteoporosis in humans and
  mice. *Nat Genet* 2019;51:258–66. — heel eBMD
