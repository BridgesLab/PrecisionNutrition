# Does serum calcium modify the LDL-C → MACE relationship?
### Progress-report summary — competing-risks interaction analysis (2026-06-08)

## Objective
Test whether the effect of LDL-C on major adverse cardiovascular events (MACE) is
modified by calcium status — the central hypothesis of the calcium–cholesterol
project. This extends the v3 main-effects analysis (in which time-averaged LDL-C
showed **no independent association** with MACE after full adjustment) by adding a
calcium × LDL-C interaction term.

## Methods
- **Design:** Fine-Gray subdistribution-hazard competing-risks model (death without
  MACE treated as the competing event), built on the v3 framework.
- **Exposure:** time-averaged LDL-C (trapezoidal AUC ÷ follow-up years), per 10 mg/dL.
- **Effect modifier:** **ionized calcium** (bioactive fraction; no albumin correction
  required; ~20,089 patients vs only ~4,009 for total/serum calcium in this cohort),
  characterized at baseline (nearest value within ±730 days of first LDL-C), per
  0.5 mg/dL.
- **Covariates:** age, sex, race/ethnicity, ever-statin, diabetes, hypertension.
- **Cohort:** patients with ≥2 LDL-C measurements and a baseline ionized-calcium value.
  **N = 367; 154 MACE events; 61 competing deaths** (events-per-variable 12.8).
- **Interaction assessed three ways** (continuous product term, calcium-tertile product,
  fully stratified models), with a Cox likelihood-ratio test and a time-averaged-calcium
  re-analysis as sensitivity checks.

## Key results

**Main effects (no interaction):** neither LDL-C (sub-HR 1.02 per 10 mg/dL, 95% CI
0.98–1.07, p=0.29) nor ionized calcium (sub-HR 1.04 per 0.5 mg/dL, p=0.59) was
independently associated with MACE — consistent with the null LDL-C main effect in v3.

**LDL-C × ionized-calcium interaction (primary, continuous):**
sub-HR ratio **1.032 per [0.5 mg/dL calcium × 10 mg/dL LDL-C]** (95% CI 1.002–1.063,
**p = 0.039**). Direction: the LDL-C effect on MACE becomes more harmful as ionized
calcium rises.

**By ionized-calcium tertile** (LDL-C sub-HR per 10 mg/dL within each tertile):

| Ionized calcium tertile | LDL-C sub-HR | 95% CI | p |
|---|---|---|---|
| T1 (low, 3.1–4.7 mg/dL) | 0.99 | 0.93–1.06 | 0.74 |
| T2 (mid, 4.8–5.1 mg/dL) | 1.03 | 0.95–1.12 | 0.45 |
| T3 (high, 5.2–7.8 mg/dL) | **1.06** | 1.00–1.13 | **0.049** |

LDL-C is unassociated with MACE at low/mid calcium and reaches nominal significance
only in the highest calcium tertile. Fully stratified models agree (T3 sub-HR 1.07,
95% CI 1.01–1.13, p=0.022).

**Cox likelihood-ratio test for the interaction:** χ²=3.88, df=1, **p=0.049**
(interaction HR 1.037, 95% CI 1.001–1.073).

## Robustness / caveats
- **Borderline and power-sensitive.** At a narrower ±365-day calcium window
  (N=219, 97 events) the interaction was stronger (sub-HR ratio 1.061, p=0.0036);
  at the better-powered ±730-day window it attenuated to p≈0.04–0.05, and the
  tertile interaction *term* itself was no longer significant (p=0.11). A robust
  effect would be expected to strengthen, not weaken, with more data — so this
  signal is fragile.
- **Not replicated with time-averaged calcium.** Using time-averaged rather than
  baseline ionized calcium, the interaction disappeared (sub-HR ratio 0.997,
  p=0.79). The signal is specific to *baseline* calcium.
- **Robust to comorbidity/severity adjustment.** Ionized calcium is preferentially
  ordered in acute/inpatient settings, raising concern that baseline calcium marks
  illness severity rather than calcium physiology. Adjusting for Charlson comorbidity
  score and renal disease left the interaction essentially unchanged:

  | Model | Calcium×LDL-C sub-HR | 95% CI | p |
  |---|---|---|---|
  | Primary (no severity adj.) | 1.032 | 1.002–1.063 | 0.039 |
  | + Charlson score + renal disease | 1.032 | 1.001–1.064 | 0.046 |
  | + Charlson × LDL-C interaction | 1.033 | 1.001–1.066 | 0.040 |

  Charlson score was itself an independent predictor of MACE (sub-HR 1.21/point,
  p=0.018) and renal disease was present in 50% of the cohort — i.e. these were
  informative covariates, yet they did not account for the calcium interaction.
  Severity did **not** modify the LDL-C effect (Charlson×LDL-C p=0.36). This argues
  against *chronic*-comorbidity confounding as the explanation.
- **Acute-acuity confounding remains untestable.** The data carry no encounter-type
  (inpatient/ED/ICU) field, so the purely acute index-state cannot be adjusted for;
  the Charlson proxy is also sparse (median 0). Acute-illness confounding at the
  calcium draw therefore cannot be fully excluded.

## Interpretation
There is a **direction-consistent signal** that LDL-C is more harmful at higher
ionized calcium, concentrated in the top calcium tertile — concordant across three
baseline-calcium specifications and a Cox LRT (all p≈0.02–0.05), and **robust to
adjustment for chronic comorbidity burden and renal disease**. It is, however, still
**statistically borderline** (p≈0.04 when adequately powered), **absent for
time-averaged calcium**, and not yet de-confounded from *acute* illness state at the
index measurement. Best characterized as a **promising, severity-robust but
not-yet-confirmatory** finding worth pursuing.

## Planned next steps
1. **Total serum calcium concordance** — re-run the same Fine-Gray interaction model
   substituting albumin-corrected total serum calcium for ionized calcium. This is
   lower-powered (~4,009 patients vs 20,089) and may well be null, but it is a natural
   sensitivity check since total calcium is the more conventional clinical measure and
   is more likely to reflect outpatient draws. (Requires first adding an albumin branch
   to `labs_cleaning.qmd`, then corrected = total + 0.8 × (4.0 − albumin).)
2. **Osteoporosis outcome** — parallel Fine-Gray competing-risks analysis flipping the
   event (osteoporosis diagnosis as primary, MACE as competing event), same interaction
   framework. This is arguably the more mechanistically motivated analysis given the
   Mendelian-randomization work.
3. **Acute-state confounding** — obtain an encounter-type (inpatient/ED/ICU) field in a
   future pull to directly test acute index-state confounding, the one alternative
   explanation not yet excluded (the Charlson proxy controls only chronic burden).
4. **Power** — pre-specify the baseline-calcium window and consider a larger
   calcium-enriched data pull to move the borderline interaction toward definitive power.
