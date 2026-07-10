# Multivariable MR: Cholesterol → BMD → Calcium Mediation

## Objective

Test whether **bone mineral density (BMD)** mediates the causal effect of
cholesterol on **serum calcium**, and whether this mediation is specific to the
**HMGCR / mevalonate** pathway rather than the **PCSK9 / LDL-receptor** or
**NPC1L1 / absorption** pathways.

**Primary hypothesis.** HMGCR cis-instruments show strong BMD-mediated effects
on calcium; PCSK9 and NPC1L1 instruments show weak or null mediation.

This document is the frozen specification for `mvmr_analyses.qmd`. It records the
decisions made at the outset so the analysis is reproducible and reviewable.

---

## Locked-in configuration (decided 2026-07-09)

| Choice | Decision | Rationale |
|---|---|---|
| **Calcium outcome** | `ebi-a-GCST90025990` (Barton 2021, UKB, n≈400,792), queried via OpenGWAS | Matches the existing `drug_target_mr_ukb.qmd` primary; rsID-based API query avoids local positional-join issues and gives full coverage of the pooled instrument set. |
| **Sample overlap** | UKB configuration is **primary**; overlap documented, not avoided | LDL-C (`ieu-b-110`), Morris heel BMD (`ebi-a-GCST006979`) and Barton calcium are all UKB → effectively one-sample MVMR. With F ≫ 10 the one-sample bias is toward the null, so a positive mediation result is conservative. Total-cholesterol (GLGC `ebi-a-GCST90025953`) is run as a non-overlapping sensitivity exposure. |
| **Mediation estimator** | **Both** difference-method MVMR (`mr_mvivw`) **and** two-step product-of-coefficients, reported side by side | Cis drug-target exposures carry few SNPs, which destabilises the difference-method *direct* effect. The two-step estimator is more stable for small cis instrument sets; agreement across the two strengthens the causal claim. |
| **MVMR robust estimator** | `MendelianRandomization::mr_mvivw` (IVW multiplicative random effects) as primary. `mv_raps` **does not exist** in `MendelianRandomization`; true multivariable MR-RAPS is `GRAPPLE::grappleRobustEst`, wired as an optional robustness hook. Univariable legs use TwoSampleMR `mr_raps`. | Honest package usage; keeps the requested RAPS robustness where it is actually implemented. |
| **Deliverable** | Commented Quarto `.qmd`; the user renders with their OpenGWAS token | Matches the rest of the repository (cached chunks, Michigan palette, `kable`). |

---

## Data sources

| Role | Trait | Source | ID | n |
|---|---|---|---|---|
| Exposure (primary) | LDL-C | UKB Neale (OpenGWAS) | `ieu-b-110` | ≈440,546 |
| Exposure (sensitivity) | Total cholesterol | GLGC meta (OpenGWAS) | `ebi-a-GCST90025953` | non-UKB |
| Mediator | Heel eBMD | Morris 2019 UKB (OpenGWAS) | `ebi-a-GCST006979` | 426,824 |
| Outcome | Serum calcium | Barton 2021 UKB (OpenGWAS) | `ebi-a-GCST90025990` | 400,792 |

**Drug-target cis windows (GRCh37, ±500 kb):** HMGCR (chr5:74,632,993–74,657,941),
PCSK9 (chr1:55,505,221–55,530,525), NPC1L1 (chr7:44,552,971–44,604,640).
Clumping: HMGCR allele score at r² < 0.30 (Swerdlow 2015); PCSK9 / NPC1L1 strict
at r² < 0.001. Cis p-value threshold 5e-8 (primary), relaxed to 1e-5 as sensitivity.

---

## Estimands & mediation algebra

For a given cholesterol exposure X, mediator M = BMD, outcome Y = calcium:

- **Total effect** `τ`: univariable MR of X → Y.
- **Direct effect** `θ`: MVMR of (X, M) → Y — coefficient on X.
- **Indirect (mediated) effect**:
  - *Difference method:* `τ − θ`.
  - *Two-step product of coefficients:* `a × b`, where
    `a` = MR(X → M) and `b` = MVMR(M | X → Y) (or univariable M → Y when X and M
    instruments are near-disjoint, as for cis targets).
- **Proportion mediated** `P_M = (τ − θ) / τ` (difference) or `a·b / τ` (product).
- **SE of the indirect effect** via the delta method:
  `SE(a·b) = sqrt(a² · SE(b)² + b² · SE(a)²)`.

Proportion-mediated is reported only when `τ` is directionally consistent and
distinguishable from zero, and is bounded-checked (values outside [0, 1] flagged,
not silently clipped).

---

## Workflow

### Part 1 — Drug-target MVMR (HMGCR, PCSK9, NPC1L1)
For each target: (1) extract cis cholesterol instruments; (2) extract genome-wide
independent BMD instruments (mediator); (3) pool rsIDs, fetch cholesterol, BMD and
calcium effects for the union, harmonise; (4) univariable cis → calcium (`τ`);
(5) difference-method MVMR (`θ`, `P_M`); (6) two-step `a·b`; (7) robustness.
Checkpoint each target to CSV.

### Part 2 — All-SNP MVMR (power / robustness)
Genome-wide LDL-C instruments (`ieu-b-110`) + BMD via TwoSampleMR
`mv_extract_exposures()` → `mv_harmonise_data()` → `mv_multiple()` and `mr_mvivw`.
Same mediation algebra. Expected to mirror HMGCR with more power.

### Part 3 — Comparison table + figures
Side-by-side table (N SNPs, total, direct, % mediated, p-direct) across all
pathways; forest plot (total vs direct by pathway); mediation bar plot (% mediated).

### Part 4 — eQTL MVMR (functional mechanism, optional)
If HMGCR eQTLs are supplied (GTEx liver / eQTLGen blood): sequential decomposition
SNP → expression → BMD → calcium (Models 1–3). Skipped gracefully if absent.

### Robustness (every model)
MR-PRESSO global + outlier test; MR-Egger intercept; Cochran's Q / I²; per-SNP and
per-exposure conditional F; leave-one-out. IVW-MRE vs (univariable) RAPS consistency.

---

## Outputs

- `results/mvmr_results_summary.csv` — Pathway, N_SNPs, Model, Effect_Type, Beta, SE, CI_lower, CI_upper, P_value, Pct_Mediated
- `results/mvmr_pleiotropy_tests.csv` — MR-PRESSO + Egger intercept per model
- `results/mvmr_heterogeneity.csv` — Cochran's Q, I², conditional F per model
- `results/mvmr_leaveoneout.csv` — LOO direct effects
- Figures: forest plot, mediation bar plot, pathway diagram (written to `figures/`)

## Interpretation

- HMGCR 70–90 % mediated through BMD → mevalonate hypothesis supported.
- PCSK9 / NPC1L1 < 20 % mediated → LDL-receptor / absorption not the main route.
- All-SNP pattern ≈ HMGCR → robust across genetic architecture.
