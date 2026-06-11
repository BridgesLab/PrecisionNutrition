# MACE/Osteoporosis LDL-C Analysis — Session Notes
# For continuation in a new chat session
# Last updated: 2026-04-06 (v3 results complete)

## Project Overview

Estimating hazard ratios for LDL-C on:
1. **Primary outcome**: MACE (major adverse cardiovascular events) — proof of principle
2. **Secondary outcome**: LDL-C effects on osteoporosis
3. **Future**: Serum calcium modification of the LDL-C → MACE relationship (interaction term)

Scientific goal is to understand the **biology/physiology** of LDL-C exposure on
cardiovascular and skeletal outcomes — not limited to demonstrating a positive
association. The analysis should be methodologically bulletproof for peer review.

---

## Data Source

**University of Michigan Precision Health DataDirect** — primary care cohort.
Full cohort: ~201,073 patients with at least one LDL-C measurement (2000–present).
Pulled in 3 batches (males / females_unmarried / females_married), combined via `combine_batches_v2.R`.

### Key files (all in `combined_data/`)
| File | Contents |
|---|---|
| `DemographicInfo.csv` | One row per patient — the **cohort base** (201,073 patients) |
| `DiagnosesCleaned.csv` | MACE and osteoporosis diagnoses only — patients absent here are controls |
| `LabResultsCleaned.csv` | All lab results — filter by `test_name == "LDL-C"` for LDL-C |
| `EncounterAll.csv` | All encounters — used to define censoring date for controls |
| `MedicationOrdersCleanedStatins.csv` | Pre-cleaned statin intervals with `period_start`, `period_end`, `intensity` |
| `ComorbiditiesOnset.csv` | Onset dates for diabetes, hypertension, obesity, renal failure |
| `MichiganDeathIndex.csv` | Death dates from Michigan Death Index (supplements EHR death dates) |

### Key column names
- Patient ID: `DeID_PatientID` (character — must coerce explicitly with `as.character()`)
- Lab date: `DeID_COLLECTION_DATE` (primary), `DeID_AdmitDate` (fallback) — format varies, parse with multiple orders
- Encounter date: `DeID_AdmitDate` in encounter file — format `mdy_hm()`
- MACE onset: `MACE.onset` in diagnosis file — format `ymd()` (date only)
- Osteoporosis onset: `Osteoporosis.onset` in diagnosis file — format `ymd()`
- LDL-C value: `value` column in lab file
- Demographic: `GenderCode` (sex), `RaceCode` (race), `EthnicityCode` (ethnicity)
  - GenderCode: `F` (n=109,410), `M` (n=91,653), `U` (n=10)
  - RaceCode: `C`=Caucasian (150,655), `AA`=African American (19,851), `A`=Asian (18,004),
    `O`=Other (7,346), `D`=Declined (1,388), `U`=Unknown (1,564), `AI`=American Indian (748),
    `P`=Pacific Islander (149), NA (1,368)
  - EthnicityCode: `NonHL`=Not Hispanic (180,630), `HL`=Hispanic/Latino (7,948),
    `U`=Unknown (5,516), `D`=Declined (1,849), NA (5,130)
  - Mapping: HL ethnicity overrides race → "Hispanic"; then C→White, AA→Black, A→Asian,
    all others→Other/Unknown
- Death: `DeID_DeceasedDate` (EHR), `DeID_MDIDeceasedDate` (Michigan Death Index)
- Statin intervals: `period_start`, `period_end`, `intensity` (low/moderate/high)

---

## Analytical Approach (v3 — Unified Competing Risks)

### Previous approaches (superseded)

The analysis went through several iterations:

1. **v1 — Time-varying spot LDL-C with stratification** (`ldlc_mace_analysis.qmd`):
   Used tmerge with LOCF and stratified Cox models (strata on statin, age category,
   period, cardiometabolic stratum). Problems: paradoxical inverse association from
   statin confounding by indication, PH violations, LOCF staleness bug for
   single-measurement patients, BMI missingness dropping 39% of sample.

2. **v2 — Cumulative LDL-years** (`ldlc_mace_cumulative.qmd`):
   Raw cumulative AUC (mg/dL × years) as exposure. Fixed PH violation but raw
   cumulative is confounded with follow-up duration (longer-lived patients
   mechanically accumulate more LDL-years). Time-averaged LDL-C (AUC ÷ years)
   was identified as the correct metric in the osteoporosis script.

3. **v3 — Unified competing risks (CURRENT)**: See below.

### Current approach: unified competing risks framework

**Rationale for the redesign (2026-04-06):**

The previous analyses had accumulated inconsistencies that made cross-outcome
comparison impossible: different exposure metrics (spot vs cumulative vs
time-averaged), different model frameworks (Cox vs Fine-Gray), different
covariate strategies (stratification vs regression), and different cohort
restrictions. The v3 redesign harmonises everything into a single framework.

**Key design decisions and rationale:**

1. **Primary exposure: time-averaged LDL-C (trapezoid AUC ÷ follow-up years), per 10 mg/dL**
   - Raw cumulative AUC is confounded with follow-up duration — patients who
     survive longer mechanically accumulate more LDL-years, creating collinearity
     with time-on-study that distorts the HR
   - Time-averaged removes this dependency: interpretation is "per 10 mg/dL higher
     average lifetime LDL-C, the hazard changes by X"
   - Trapezoid rule for AUC is primary; step/LOCF is sensitivity
   - Biologically motivated: captures sustained exposure rather than any single point

2. **Fine-Gray competing risks (primary model)**
   - MACE: death without MACE as competing event (two-state)
   - Osteoporosis: MACE and death as competing events (three-state)
   - Standard Cox treats death as non-informative censoring, which biases the
     estimate when LDL-C affects both the outcome and mortality risk
   - Fine-Gray models the subdistribution hazard, accounting for competing events
   - Standard Cox included as sensitivity for comparison

3. **Covariates as regression terms (not stratification)**
   - Previous MACE model stratified on statin, age category, period, and
     cardiometabolic stratum — this prevented estimating the statin effect,
     consumed degrees of freedom in sparse strata (some had only 13–19 events),
     and was incompatible with `cmprsk::crr()` which doesn't support stratification
   - Moving covariates into the regression enables direct HR comparison across
     model specifications and across outcomes
   - Common covariate set for both outcomes: age at baseline, sex, race/ethnicity,
     ever/never statin, diabetes, hypertension

4. **Ever/never statin (not time-varying)**
   - People generally go on statins when indicated and stay on them
   - For the primary Fine-Gray model, ever/never is a fixed covariate —
     `crr()` does not support time-varying covariates
   - The time-varying Cox sensitivity uses on/off statin at each interval,
     providing a check on whether the simpler ever/never is adequate

5. **Race/ethnicity included**
   - Combined from `RaceCode` + `EthnicityCode` using NIH convention:
     Hispanic/Latino ethnicity overrides race category
   - Categories: White (reference), Black, Asian, Hispanic, Other/Unknown
   - GenderCode needs verification — current regex maps `^f` → Female, `^m` → Male;
     check actual values if the mapping produces zero-variance columns

6. **≥2 LDL-C measurements (primary cohort)**
   - Single-measurement patients give a point estimate masquerading as a
     time-average — the trapezoid AUC for a single point is just
     `LDL × follow-up`, making the time-average equal to the single measurement
   - ≥2 measurements primary; all patients as sensitivity

7. **Death date from two sources**
   - EHR death (`DeID_DeceasedDate`) and Michigan Death Index (`DeID_MDIDeceasedDate`)
   - Unified as the earliest of the two
   - MDI adds ~2,400 additional deaths not in EHR
   - Used only for competing event coding, not for extending follow-up
     (censoring remains at last encounter to avoid informative censoring)

### Model specifications

**Primary model (both outcomes):**
Fine-Gray subdistribution hazard with time-averaged LDL-C (trapezoid, per 10 mg/dL),
adjusted for age at baseline, sex, race/ethnicity, ever/never statin, diabetes,
hypertension. Cohort restricted to ≥2 LDL-C measurements.

**Sensitivity analyses (run in each script):**

| # | Analysis | Purpose |
|---|---|---|
| 1 | Standard Cox (same covariates) | Quantify impact of competing risks adjustment |
| 2 | Baseline LDL-C (Fine-Gray) | Compare instantaneous vs cumulative exposure |
| 3 | Time-varying LOCF via tmerge (Cox) | Most granular time-updating; uses time-varying statin |
| 4 | Step/LOCF AUC (osteo only) | Check sensitivity to AUC calculation method |
| 5 | Statin-naive subgroup (Fine-Gray) | Remove statin confounding entirely |
| 6 | All patients incl. single measurement | Check if ≥2 restriction changes results |
| 7 | Two-state vs three-state (osteo only) | Quantify contribution of death as competing event |

### R packages required
```r
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(cmprsk)         # Fine-Gray competing risks
library(broom)
library(knitr)
```

---

## Cohort Definition

```
Full demographic cohort:                 201,073
No LDL-C data:                         ~177,667
MACE/Osteo=TRUE with missing onset:      ~2,867
Negative or zero follow-up:               ~556
──────────────────────────────────────────────
Final analytic cohort:                   ~20,015
  ≥2 LDL-C measurements:                ~5,686  (primary for osteoporosis)
```

**Time zero**: date of first LDL-C measurement
**Censoring**: last encounter date (from EncounterAll.csv)
**Events**: MACE onset / Osteoporosis onset
**Competing events**: death (both outcomes); MACE (osteoporosis only)

Note: exact numbers differ slightly between MACE and osteoporosis cohorts because
different exclusion criteria apply (MACE-missing-onset vs osteo-missing-onset).

---

## Code Files

| File | Status | Purpose |
|---|---|---|
| `combine_batches_v2.R` | Active | Combines 3 data pull batches into `combined_data/` |
| `ldlc_mace_competing_risks.qmd` | **ACTIVE — v3** | MACE primary analysis with competing risks |
| `ldlc_osteoporosis_competing_risks.qmd` | **ACTIVE — v3** | Osteoporosis primary analysis with competing risks |
| `ldlc_calcium_interaction_mace.qmd` | **NEW — feasibility-limited** | Serum calcium × LDL-C effect modification on MACE (Fine-Gray), built on v3 framework. Calcium is sparse in the single combined_data pull; see the in-script "Calcium Feasibility Diagnostic" before trusting estimates. |
| `ldlc_mace_analysis.qmd` | Superseded (v1) | Time-varying spot LDL-C with stratification |
| `ldlc_mace_baseline.qmd` | Superseded (v1) | Fixed baseline LDL-C |
| `ldlc_mace_cumulative.qmd` | Superseded (v2) | Raw cumulative LDL-years |
| `ldlc_osteoporosis_cumulative.qmd` | Superseded (v2) | Osteoporosis with multiple exposure metrics |
| `step1_data_qc.R` | Superseded | Original QC script (surgical cohort) |
| `step2_tmerge_dataset_v2.R` | Superseded | Original tmerge script (surgical cohort) |
| `step3_cox_model.R` | Superseded | Original Cox model script (surgical cohort) |

---

## Current Results (v3 — Fully Adjusted Competing Risks)

Results from `ldlc_mace_competing_risks.qmd` and `ldlc_osteoporosis_competing_risks.qmd`,
run 2026-04-06. Primary model: Fine-Gray, time-averaged LDL-C per 10 mg/dL,
adjusted for age, sex, race/ethnicity, ever statin, diabetes, hypertension.
Cohort: ≥2 LDL-C measurements.

### MACE — Primary Model (Fine-Gray, death competing)

| Term | Sub-HR | 95% CI | p |
|---|---|---|---|
| **Time-avg LDL-C (per 10 mg/dL)** | **1.004** | **0.985–1.023** | **0.67** |
| Age at baseline | 1.032 | 1.027–1.037 | <0.001 |
| Female (vs Male) | 0.885 | 0.773–1.012 | 0.075 |
| Black (vs White) | 1.064 | 0.853–1.326 | 0.58 |
| Asian (vs White) | 0.519 | 0.375–0.718 | <0.001 |
| Hispanic (vs White) | 0.806 | 0.540–1.204 | 0.29 |
| Ever statin | 0.746 | 0.637–0.874 | <0.001 |
| Diabetes | 1.220 | 1.060–1.403 | 0.005 |
| Hypertension | 1.888 | 1.512–2.357 | <0.001 |

**N = 5,668 patients; 1,032 MACE events.**

### MACE — All Model Specifications

| Model | HR | 95% CI | p |
|---|---|---|---|
| Fine-Gray — time-avg (PRIMARY) | 1.004 | 0.985–1.023 | 0.67 |
| Fine-Gray — baseline LDL-C | 1.001 | 0.984–1.018 | 0.92 |
| Fine-Gray — statin-naive | 0.972 | 0.936–1.011 | 0.15 |
| Fine-Gray — all patients | 0.987 | 0.975–0.998 | 0.024 |
| Cox — time-avg LDL-C | 1.003 | 0.985–1.021 | 0.74 |
| Cox — baseline LDL-C | 1.000 | 0.984–1.016 | 1.00 |
| Cox — time-varying LOCF 1825d | 1.002 | 0.982–1.022 | 0.85 |

### Osteoporosis — Primary Model (Fine-Gray, MACE + death competing)

| Term | Sub-HR | 95% CI | p |
|---|---|---|---|
| **Time-avg LDL-C (per 10 mg/dL)** | **0.985** | **0.962–1.009** | **0.21** |
| Age at baseline | 1.055 | 1.048–1.061 | <0.001 |
| Female (vs Male) | 6.035 | 4.941–7.372 | <0.001 |
| Black (vs White) | 0.662 | 0.461–0.951 | 0.026 |
| Asian (vs White) | 0.873 | 0.619–1.232 | 0.44 |
| Hispanic (vs White) | 0.584 | 0.317–1.076 | 0.084 |
| Ever statin | 0.646 | 0.530–0.787 | <0.001 |
| Diabetes | 0.782 | 0.650–0.941 | 0.009 |
| Hypertension | 0.992 | 0.762–1.292 | 0.95 |

**N = 5,784 patients; 595 osteoporosis events.**

### Osteoporosis — All Model Specifications

| Model | HR | 95% CI | p |
|---|---|---|---|
| Fine-Gray 3-state — time-avg (PRIMARY) | 0.985 | 0.962–1.009 | 0.21 |
| Fine-Gray 3-state — baseline LDL-C | 0.989 | 0.969–1.010 | 0.29 |
| Fine-Gray 3-state — step/LOCF | 0.987 | 0.964–1.011 | 0.28 |
| Fine-Gray 3-state — statin-naive | 0.999 | 0.954–1.046 | 0.96 |
| Fine-Gray 3-state — all patients | 0.992 | 0.980–1.004 | 0.19 |
| Fine-Gray 2-state — MACE competing only | 0.983 | 0.960–1.007 | 0.16 |
| Cox — time-avg LDL-C | 0.981 | 0.959–1.004 | 0.10 |
| Cox — baseline LDL-C | 0.987 | 0.967–1.007 | 0.19 |
| Cox — time-varying LOCF 1825d | 0.984 | 0.957–1.011 | 0.24 |

### Interpretation

**LDL-C shows no independent association with either MACE or osteoporosis**
after full covariate adjustment. This null finding is consistent across all
model specifications (7 for MACE, 9 for osteoporosis), exposure parameterisations
(time-averaged, baseline, time-varying), and model frameworks (Fine-Gray, Cox).

**Key covariate effects validate the model:**
- **MACE**: Hypertension (Sub-HR 1.89), diabetes (1.22), and statin use (0.75)
  are all strongly associated in expected directions. Age, sex, and race effects
  are biologically plausible. This rules out model misspecification as an
  explanation for the null LDL-C finding.
- **Osteoporosis**: Female sex (Sub-HR 6.04), age, and Black race (0.66, protective)
  are all well-established. Statin use is protective (0.65), consistent with
  emerging evidence for statin bone effects. Diabetes is protective (0.78),
  consistent with higher BMD in type 2 diabetes.

**Why does the v2 osteoporosis finding (Sub-HR 1.05, p<0.001) disappear?**
The v2 model adjusted only for age and sex. Adding statin use and diabetes
completely absorbed the apparent LDL-C → osteoporosis association:
- Statin users had both lower LDL-C AND lower osteoporosis risk (Sub-HR 0.65),
  creating a spurious positive LDL-C–osteoporosis correlation
- Diabetic patients had higher LDL-C AND lower osteoporosis diagnosis rates
  (Sub-HR 0.78), reinforcing the same confounding direction
- This demonstrates the critical importance of adequate confounder adjustment
  and is itself a publishable methodological finding

**The statin paradox in the MACE analysis:**
The statin-naive subgroup (Sub-HR 0.97, p=0.15) trends slightly protective —
the opposite of what a naive "remove confounding" approach would predict. This
is consistent with collider/selection bias (Berkson's paradox): conditioning on
statin-naive selects for patients where high LDL-C and high MACE risk do not
co-occur, creating a spurious inverse association.

**Power considerations:**
The primary MACE CI (0.985–1.023 per 10 mg/dL) can rule out effects larger than
~2.3% per 10 mg/dL. The expected causal effect from CTT meta-analyses (RCT data)
is approximately 5–10% per 10 mg/dL. The null finding likely reflects:
(a) the dominant statin confounding pathway operating through non-LDL mechanisms,
(b) limited power with ~1,032 events in the ≥2 measurement cohort, and
(c) the observational setting where treatment-by-indication obscures the causal effect.

---

## Previous Results (from v1/v2 — for reference only)

These results are from earlier model specifications and should NOT be cited.
They are retained here to track the evolution of the analysis.

### v1: Time-varying spot LDL-C (per 1 mg/dL, stratified Cox)
| Model | HR | 95% CI | p |
|---|---|---|---|
| Unadjusted (730d) | 0.9966 | 0.9955–0.9977 | <0.001 |
| Age+BMI adjusted (730d) | 0.9991 | 0.9976–1.0007 | 0.27 |

Paradoxical inverse association driven by statin confounding by indication.

### v2: Cumulative LDL-years (per 40 mg/dL-years, stratified Cox)
| Model | HR | 95% CI | p |
|---|---|---|---|
| Unadjusted (730d) | 1.0112 | 1.008–1.014 | <0.001 |
| Stratified + statin (730d) | 1.0038 | 0.998–1.010 | 0.24 |

### v2: Osteoporosis Fine-Gray (per 10 mg/dL, age+sex adjusted)
| Model | Sub-HR | 95% CI | p |
|---|---|---|---|
| Fine-Gray baseline LDL-C | 1.0359 | 1.013–1.059 | 0.002 |
| Fine-Gray time-avg trapezoid | 1.0486 | 1.024–1.074 | <0.001 |
| Cox time-avg trapezoid | 1.0000 | 0.976–1.024 | 0.997 |

Notable: Fine-Gray showed significant positive association for osteoporosis that
was completely absent in standard Cox — demonstrates the importance of competing
risks adjustment.

---

## Known Issues and TODOs

### Near-term

- [~] **Serum calcium interaction**: SCRIPT BUILT (`ldlc_calcium_interaction_mace.qmd`,
  2026-06-08). Modifier = baseline serum calcium (`test_name == "Serum Calcium"`,
  mg/dL; nearest to t0 within ±365d); time-averaged calcium as sensitivity. Three
  views: continuous LDL×Ca product (primary Wald test), LDL×Ca-tertile (sub-HR per
  tertile via linear combos), fully stratified Fine-Gray. Cox LRT + time-avg-Ca as
  sensitivities.
  - ⚠ **COHORT LESSON (2026-06-08)**: only 66/20,015 base-cohort patients had a
    baseline calcium, so the analytic cohort collapsed to **N=19, 3 MACE events**.
    The resulting "interaction p=0.0015, sub-HR 65.7" was a QUASI-COMPLETE-
    SEPARATION ARTIFACT (sub-HRs of 0 and 1e+17, Cox LRT p=1 with Inf CIs) — NOT
    signal. Discarded.
  - ⚠ **DATA REALITY (2026-06-08)**: there is only ONE combined_data (the
    2026-03-23 pull). The data-pull notes describing a separate 50,225 "Calcium
    and Cholesterol" MGI extraction appear to be a documentation error — no
    distinct calcium-rich cohort exists to switch to. Serum calcium (~4,000
    patients with any value) is the binding constraint.
  - ✅ **COVERAGE FIX (2026-06-08)**: `ca-chol-correlations.qmd` shows IONIZED
    calcium has ~20,089 patients vs serum/total ~4,009 (LDL-C ~23,406). Script
    now has `CALCIUM_TEST` config defaulting to **"Ionized Calcium"** (bioactive,
    no albumin correction, ~5× coverage; but more inpatient/acute context).
    Serum calcium = sensitivity. Scales differ (ionized ~4.5–5.3, total ~8.5–10.5
    mg/dL) so guardrails + CAL_SCALE auto-switch by measure. The project's own
    correlation work (Gamma/brms) already used ionized calcium.
  - **CURRENT STATE**: script has `CALCIUM_TEST`/`DATA_DIR` config, an events-per-
    variable guard, and a "Calcium Feasibility Diagnostic" table (N / MACE events
    / EPV crossing LDL requirement × baseline-calcium window).
  - 🔵 **FIRST REAL RESULT (2026-06-08, ionized calcium, ≥2 LDL, ±365d window,
    N=219/97 MACE events)**: a SIGNIFICANT positive LDL-C × baseline-ionized-Ca
    interaction on MACE. Continuous interaction sub-HR ratio 1.061 (1.020–1.105)
    per [0.5 mg/dL Ca × 10 mg/dL LDL], p=0.0036. Monotonic by tertile: LDL sub-HR
    T1(low) 0.96 ns → T2 1.07 ns → T3(high) 1.09–1.10, p≈0.001. Cox LRT χ²=7.47,
    p=0.0063. LDL main effect at mean Ca null (~1.03), consistent with v3.
    Interpretation: LDL-C harmful only at high ionized calcium.
  - ⚠ **CAVEATS (do not over-claim)**: (1) Time-averaged-Ca sensitivity does NOT
    replicate (ldl_x_cavg 0.986, p=0.31, N=176) — interaction is specific to
    BASELINE calcium. (2) Ionized calcium is ordered in acute/inpatient settings,
    so baseline ionized Ca may mark acuity at index → possible confounding by
    illness severity, not stable biology. (3) EPV was 8.1 (<10) at ±365d.
  - 🔵 **POWERED RESULT (730d, N=367/154ev/EPV 12.8)**: interaction attenuated vs
    365d but persists — ldl_x_cac 1.032 (1.002–1.063), p=0.039; tertile T3 LDL
    sub-HR 1.06 (1.00–1.13) p=0.049; Cox LRT p=0.049. Borderline.
  - ✅ **SEVERITY ADJUSTMENT (Sensitivity C)**: added Charlson baseline score +
    renal disease (no encounter-type field exists for true acute acuity). Interaction
    UNCHANGED across Primary→C1→C2 (1.0317→1.0319→1.0330, p~0.04). Charlson itself
    predicted MACE (1.21/pt, p=0.018) and renal disease was 50% prevalent, so the
    covariates were live — argues AGAINST chronic-severity confounding. Charlson×LDL
    null (p=0.36). Caveat: Charlson baseline sparse (median 0); acute index-state
    still untestable without encounter-type data.
  - **STATUS**: promising, severity-robust, but NOT confirmatory — still borderline
    (p~0.04) and absent for time-averaged calcium. NEXT: serum-calcium concordance
    check (+albumin correction); pursue encounter-type field for acute-acuity test;
    pre-specify window; consider larger calcium-enriched pull for definitive power.
  - TODO: from the diagnostic, pick a defensible (window, LDL) combo or document
    underpowering; consider albumin correction (albumin not yet in
    `labs_cleaning.qmd` test_name — would need a new branch).
- [ ] **LDL-C × diabetes and LDL-C × hypertension interactions**: the strong
  covariate effects for diabetes (MACE Sub-HR 1.22) and hypertension (MACE
  Sub-HR 1.89) suggest potential effect modification. Diabetic dyslipidemia
  produces small dense LDL (more atherogenic per unit LDL-C); hypertension
  increases endothelial shear stress amplifying LDL infiltration. A null
  LDL-C main effect could mask real effects concentrated in comorbid subgroups.
  Add interaction terms to Cox models (Fine-Gray interactions trickier but doable).
- [ ] **LDL-C × sex interaction for osteoporosis**: with female Sub-HR 6.04
  dominating the model, the LDL-C effect may differ substantially by sex.
  Test as interaction term.
- [ ] **BMI handling**: BMI was dropped from the v3 primary models due to massive
  missingness (~39% sample loss). Options: multiple imputation via `mice`,
  check if height/weight columns exist separately, or keep BMI out of
  primary and include as sensitivity in complete-case subset.
- [ ] **MACE definition review**: current ICD codes may capture treated prevalent
  disease rather than incident hard events. Consider tightening to acute MI
  (I21.x, I22.x), stroke (I63.x), CV death (I46.x, I51.x) only.
- [ ] **Osteoporosis-specific confounders**: bisphosphonates, calcium/vitamin D
  supplementation, corticosteroids, hormone therapy. These were identified as
  important but not yet incorporated. Should be added as sensitivity analyses
  to the osteoporosis script.
- [ ] **Progressive adjustment table for manuscript**: show how the osteoporosis
  LDL-C Sub-HR changes from ~1.05 (age+sex only) → ~1.00 (fully adjusted)
  as confounders are added sequentially. This demonstrates confounding and is
  a key methodological contribution.

### Longer-term / manuscript preparation

- [ ] **Formal HR comparison across outcomes**: consider interaction test
  (outcome × exposure) to formally test whether LDL-C effect differs between
  MACE and osteoporosis. Alternatively, visual comparison via side-by-side
  forest plots on the same scale.
- [ ] **Dose-response / non-linearity**: consider restricted cubic splines for
  time-averaged LDL-C to check for non-linear exposure-response. Fine-Gray
  with splines may require `FGR()` from the `riskRegression` package.
- [ ] **Statin–osteoporosis association**: the protective Sub-HR 0.65 for statins
  on osteoporosis is itself interesting. Could warrant follow-up analysis —
  is this causal (statin promotes osteoblast differentiation) or residual
  confounding (healthy-user bias)?
- [ ] **Multiple testing**: with many sensitivity analyses, consider whether
  formal multiplicity adjustment is needed or whether the sensitivity framework
  speaks for itself.
- [ ] **SES covariates**: check if DataDirect provides any socioeconomic status
  proxies (insurance type, zip code-level indices).

---

## Resolved Issues (from previous sessions)

- ✅ **Statin confounding**: resolved by including ever/never statin as covariate
  (and time-varying statin in Cox sensitivity)
- ✅ **LOCF staleness bug**: sidestepped by switching primary exposure from
  time-varying spot LDL-C to time-averaged LDL-C (no LOCF needed for primary model)
- ✅ **PH violation**: resolved in the osteoporosis Fine-Gray models; to be
  verified in MACE v3 (PH check included in scripts)
- ✅ **Cumulative LDL confounded with time**: resolved by using time-averaged
  (AUC ÷ years) instead of raw cumulative
- ✅ **Cross-outcome comparability**: resolved by harmonising exposure metric,
  covariate set, and cohort restriction across both outcomes
- ✅ **Death as competing event**: incorporated via Michigan Death Index +
  EHR death dates into Fine-Gray framework
- ✅ **Race/ethnicity**: added to models using RaceCode + EthnicityCode
- ✅ **Demographic code mapping**: DataDirect uses abbreviated codes (F/M, C/AA/A, HL/NonHL);
  fixed from regex matching (which failed on abbreviations) to exact string matching
- ✅ **v3 results reviewed (2026-04-06)**: LDL-C null for both MACE and osteoporosis
  after full adjustment. v2 osteoporosis positive finding explained by statin and
  diabetes confounding. Covariate effects validate model specification.

---

## Technical Notes

- Encounter file has 9M+ rows — loads slowly. Consider saving `last_encounter` as RDS.
- ID column `DeID_PatientID` must be coerced to `as.character()` early to avoid
  silent join failures (was a major bug source in early pipeline).
- `cmprsk::crr()` requires a numeric design matrix — factors must be dummy-coded
  manually. A `safe_crr()` wrapper is used in both scripts to automatically drop
  zero-variance columns before fitting, preventing singular matrix errors.
- **Do NOT use `with(df, crr(...))` syntax** — caused persistent scoping bugs where
  the cleaned design matrix was not visible inside `with()`. All `crr()` calls must
  pass vectors directly: `crr(ftime = df$col, fstatus = df$col, cov1 = mat, ...)`.
  Design matrices are built with `as.matrix(df %>% select(...))` not `with(df, cbind(...))`.
- Ever/never statin and ever/never diabetes/hypertension use information from the
  entire follow-up period (diagnosed at any point before event/censoring). This is
  standard for observational competing risks analyses but is not strictly a baseline
  covariate. The time-varying Cox sensitivity addresses this for statins.
- The two scripts (`ldlc_mace_competing_risks.qmd` and
  `ldlc_osteoporosis_competing_risks.qmd`) are fully self-contained and do not
  depend on each other or on superseded scripts.
