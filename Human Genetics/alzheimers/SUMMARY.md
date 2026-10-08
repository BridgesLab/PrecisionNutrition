# Body composition, glycemia and Alzheimer's disease — what we know (2026-10-08)

Status summary of the robust univariable MR pipeline (`_targets.R`, `pipeline/`) and the planned
mediation analysis (`mr_mediation.qmd`). Numbers are from `results/pipeline/` (pipeline run of
2026-10-07, commit `071fa52`); the browsable version with plots is `robust_mr_report.qmd`.
Design details and deviations are in [`pipeline/README.md`](pipeline/README.md).

## The question

Does body composition affect AD risk, and does part of lean mass's apparent protective effect run
through lower post-load (2 h) glucose? Motivating triad: **lean mass → 2 h glucose → AD**.

## Data

| Role | Dataset | Notes |
|---|---|---|
| Body composition | Appendicular lean mass (Pei 2020, `ebi-a-GCST90000025`) | **adjusted for appendicular fat mass** (collider risk) |
| | Whole-body fat-free mass (`ukb-b-13354`), fat mass (`ukb-b-19393`), BMI (`ukb-b-2303`) | UKB, MRC-IEU, SD units, not adjusted for adiposity |
| Glycemic | 2 h glucose (Chen 2021, `ebi-a-GCST90002227`) | **BMI-adjusted**. Every MAGIC 2hGlu release is (checked MAGIC's downloads page) |
| | Fasting glucose: Dupuis 2010 (`ieu-b-4761`, primary, unadjusted) and Manning 2012 (`ebi-a-GCST005186`, adjustment unclear) | |
| | T2D (Xue 2018, `ebi-a-GCST006867`) | includes UKB |
| AD | Bellenguez 2022 (`ebi-a-GCST90027158`), proxy + clinical | lacks rs429358 / rs7412 (APOE ε) |
| | Kunkle 2019 IGAP (`ieu-b-2`), clinical only, no UKB | the overlap-free, proxy-free check |
| Survival | Parental lifespan (`ebi-a-GCST006697`), Martingale residual, + = shorter life | selection axis for SlopeHunter |

## Methods, briefly

- **Main:** MR-APSS (models correlated pleiotropy and sample overlap), instruments p < 5×10⁻⁵, plus
  a **p < 5×10⁻⁸ threshold sensitivity**.
- **Sensitivity:** CAUSE (credible interval for γ is the decision rule).
- **Category-2 comparison:** MRBEE, MR-RAPS (assume InSIDE; flagged unstable only on process grounds).
- **Extra sensitivity:** MRAID, only for overlap-free pairs with ≤ 1,000 candidates. Its fits for
  lean mass / fat-free mass (3.5–5.3k candidates, 1000G LD) were **unreliable** (effects > 1 SD/SD,
  opposite signs for two FG GWAS).
- **Classical suite** with Q, I², I²_GX, F, Steiger; scatter / funnel / leave-one-out plots per pair.
- **Survival collider:** Bellenguez outcomes also run SlopeHunter-adjusted (b_SH = −0.944).
- **Latent heritable confounding:** CAUSE q·η and the MR-APSS background slope, every pair in both
  directions.
- 80 associations: exposure → mediator, adiposity → glycemic, → AD, → lifespan, and their reverses.

Scales: MR-APSS, CAUSE and MRAID are SD(outcome) per SD(exposure); classical methods and MRBEE are
per unit of each GWAS's beta (log-OR for AD). Compare across those families on sign and
significance, not magnitude.

## Findings

### 1. Lean / fat-free mass is protective for AD, robustly

| → AD (MR-APSS, SD per SD) | Bellenguez raw | Bellenguez SlopeHunter | Kunkle (clinical) |
|---|---|---|---|
| **Lean mass** | −0.017, p = 0.002 | −0.019, p = 6×10⁻⁴ | −0.041, p = 0.002 |
| **Fat-free mass** | −0.030, p = 3×10⁻⁵ | −0.022, p = 0.002 | −0.047, p = 0.005 |

- CAUSE agrees for lean mass in all three, with every 95% credible interval excluding 0
  (e.g. Kunkle γ = −0.024 [−0.048, −0.001]). For fat-free mass the CAUSE interval just reaches 0
  after adjustment and in Kunkle.
- p < 5×10⁻⁸ instruments give the same answer.
- The effect **survives survival-bias correction** (lean mass → lifespan is null, so the correction
  barely moves it) **and replicates in clinical-only, UKB-free Kunkle**.
- Leave-one-out is flat: no single SNP moves lean mass → Bellenguez by more than 0.006.
- Caveat: the Pei lean-mass GWAS is fat-mass-adjusted, but the unadjusted fat-free mass agrees.

### 2. The apparent protection from adiposity is survival bias

| → AD (MR-APSS) | Bellenguez raw | Bellenguez SlopeHunter | Kunkle |
|---|---|---|---|
| BMI | −0.023, p = 0.016 | +0.007, p = 0.48 | −0.004, p = 0.85 |
| Fat mass | −0.029, p = 0.002 | −0.002, p = 0.80 | −0.021, p = 0.35 |

Both are abolished by SlopeHunter (CAUSE also goes to ~0) and absent in Kunkle, as Step 2 found for
BMI with IVW. AD liability also predicts **lower** fat mass (reverse direction: CAUSE −0.079, CrI
excludes 0; MR-APSS −0.11), consistent with prodromal weight loss.

### 3. 2 h glucose → AD: suggestive, now supported by the main method on strong instruments

| Bellenguez | MR-APSS (5×10⁻⁵) | MR-APSS (5×10⁻⁸) | MRAID | CAUSE | IVW-MRE |
|---|---|---|---|---|---|
| raw | 0.013, p = 0.58 | **0.049, p = 0.032** | 0.030, p = 0.010 | 0.021, P(γ>0) = 0.93 | 0.10 log-OR/mmol/L, p = 0.035 |
| SlopeHunter | 0.028, p = 0.25 | **0.066, p = 0.008** | 0.044, p = 2×10⁻⁴ | 0.023, P(γ>0) = 0.94 | 0.13, p = 0.002 |

- With genome-wide-significant instruments, MR-APSS agrees with the classical methods, MRBEE,
  MR-RAPS and MRAID. The earlier null was **dilution by weak p < 5×10⁻⁵ instruments**
  (the 2hGlu GWAS has N = 63k), not the robust model discounting the signal.
- SlopeHunter adjustment strengthens every method, and the adjusted arm survives FDR.
- Kunkle is positive at every threshold but underpowered (e.g. MR-APSS 5×10⁻⁸: 0.059, p = 0.28).
- **Still conditional on the BMI adjustment of the 2hGlu GWAS** (see 5).

### 4. Fasting glucose and T2D → AD: no consistent effect

Fasting glucose (Dupuis or Manning) and T2D are null in Bellenguez under every main/sensitivity
method. The isolated Kunkle signals seen earlier (FG Dupuis by MRAID, T2D by MRAID) are not
supported by MR-APSS or CAUSE.

### 5. 2 h glucose estimates are distorted by its GWAS's BMI adjustment

| → | 2 h glucose (BMI-adjusted) | Fasting glucose (Dupuis) | T2D |
|---|---|---|---|
| BMI | **−0.092** | +0.078 | +0.21 |
| Fat mass | **−0.12** | +0.080 | +0.16 |
| Fat-free mass | **−0.14** | +0.009 | +0.042 |
| Lean mass | **−0.078** | +0.021 (ns) | +0.001 (ns) |

(MR-APSS.) Everything that raises body size gets a spurious **negative** effect on BMI-adjusted
2hGlu, the textbook collider pattern, while the same exposures raise or don't affect fasting glucose
and T2D. BMI → T2D (IVW OR ≈ 3 per SD) is the positive control and works.

Mechanism: adjusting the GWAS for BMI subtracts c × (each variant's BMI effect) from its 2hGlu
association, where c is the *observed* 2hGlu-on-BMI slope in the cohorts. The MR estimate for a
body-size exposure is then roughly (its true causal effect on 2hGlu) − c. That is zero only if the
observed slope is entirely causal; any confounded or reverse component makes c larger and the
estimate negative. So this is over-adjustment from conditioning on a heritable covariate, not
incomplete adjustment. The Dupuis fasting glucose is **not** BMI-adjusted (primary analysis), so
its positive association with BMI is the expected biology.

### 6. Latent residual heritable confounding

- **CAUSE finds no detectable shared-factor contribution anywhere:** q·η ≈ 0 with narrow intervals
  for all 72 pairs; q = 0.02–0.07, never near the 0.5 flag. Those q values sit near CAUSE's prior,
  so read this as "not detected" rather than "absent". The fraction q·η/(γ + q·η) is withheld for 45
  pairs whose denominator interval spans 0.
- **The MR-APSS background slope tracks the causal estimate** (e.g. lean mass → 2hGlu −0.078 vs
  causal −0.078; fat-free mass → Bellenguez −0.031 vs −0.030). Ω includes the causal path, so this
  is the expected upper-bound behaviour, not separable confounding. Reverse-direction slopes are
  large mostly because they scale with √(h²_outcome / h²_exposure).
- Sample-structure (overlap) slopes C_xy/C_xx (β scale) are small for most pairs but reach
  0.14–0.25 for T2D ↔ BMI / fat mass / fat-free mass, where both GWAS include UK Biobank. That's the
  overlap MR-APSS and CAUSE correct for; overlap-naive methods (IVW, Egger, RAPS…) should be
  discounted for those pairs (`overlap_C12`, `method_handles_overlap` in `results.csv`).

### 7. Survival-collider correction is valid for Bellenguez

b_SH = −0.944 is identified (raw lifespan–AD r = −0.18 on 2,243 clumped SNPs), reproduced by refit
(−0.935, SE 0.080), and lies outside its permutation null (−0.41 to +0.22, 8 permutations). Kunkle's
b_SH is not identified and is never applied.

## Where this leaves the mediation hypothesis

| Leg | Evidence |
|---|---|
| X → Y (lean mass → AD) | **Robust**: main and sensitivity methods, survival-adjusted, replicated in Kunkle |
| M → Y (2hGlu → AD) | **Suggestive**: significant with strong instruments under the main method; Kunkle underpowered |
| X → M (lean mass → 2hGlu) | **Not interpretable yet**: strong and negative, but in the direction the BMI-adjustment artefact predicts; no effect on unadjusted fasting glucose |

Mediation through 2 h glucose is neither supported nor ruled out. The decisive missing piece is a
BMI-**unadjusted** 2hGlu GWAS. With it, BMI → 2hGlu should turn positive (a built-in check the
adjustment was removed), and lean mass → 2hGlu either holds or disappears.

## Open items

1. **Unadjusted 2hGlu:** ask MAGIC directly; check Chen 2021's supplementary note on collider bias
   for lead-variant unadjusted effects; fall back to an approximate reversal of the adjustment as a
   sensitivity analysis.
2. **Mediation (`mr_mediation.qmd`):** swap in the final exposure/mediator choice once (1) is settled;
   it reuses this pipeline's instrument and SlopeHunter conventions.
3. **Optional:** more SlopeHunter permutations (8 now); a common-scale conversion of MR-APSS / CAUSE
   to per-unit estimates; MVMR with an explicit confounder, re-running CAUSE / MR-APSS to see how much
   the latent terms shrink.
4. **Housekeeping:** duplicate `apss_C12` column in `results.csv` (cosmetic; fix at table assembly).

## Bugs and pitfalls found along the way (all fixed)

- TwoSampleMR 0.7.5 `mr_ivw_mre()` doesn't floor the residual SE at 1 (SE below fixed effects when
  instruments are under-dispersed). Step 2 isn't affected.
- Kunkle APOE p-values (10⁻⁸⁸¹) underflow; now floored, not dropped.
- `tar_map()` rewrites symbols matching target names; MR-APSS calls `mvtnorm` without declaring it;
  Great Lakes' plink2 modules lack `--clump`; MR-APSS discards its LDSC jackknife delete-one values
  (recovered for the background-slope interval).
