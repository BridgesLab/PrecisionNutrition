# ANALYSIS.md — Cholesterol → Serum Calcium: mechanism tracking

**Living document.** Records what we have tested, what we found, and what remains
open, so the reasoning does not evaporate into scattered `.html.md` outputs.
Update the relevant section (and the changelog at the bottom) as new analyses land.

_Last updated: 2026-07-11._

---

## The question

Cholesterol raises serum calcium — but **how**? The original hypothesis was a
serial chain, specific to the HMGCR/mevalonate pathway:

> cholesterol ↑ → bone demineralization (BMD ↓) → calcium released into serum (↑)

The work below tests each link, and progressively **rules the bone chain out** in
favour of a bone-independent, calcium-specific mechanism.

---

## Headline conclusions (current state)

1. **Cholesterol → serum calcium is causal and directionally one-way** (forward,
   not reverse).
2. **The effect is HMGCR/mevalonate-specific**, not generic LDL-lowering — the
   HMGCR cis effect (~0.185) is ~3–5× the all-cholesterol effect (~0.05), and
   PCSK9 is informatively null.
3. **It is NOT mediated by bone.** Cholesterol lowers BMD, but BMD has no causal
   effect on serum calcium — at any skeletal site. Bone loss and calcium rise are
   **parallel** consequences of HMGCR, not a chain.
4. **It is NOT an albumin/assay artefact**, and **NOT** the PTH / vitamin-D /
   intestinal-absorption / bone-resorption axes (phosphate is unmoved; vitamin D
   moves the *wrong* way).
5. **Leading hypothesis: a calcium-specific, plausibly renal (tubular) mechanism.**
   Not yet tested directly.
6. **PTH remains genuinely untested** (no usable genetic instrument found yet) —
   a gap, not an exclusion.

---

## Evidence chain

### 1. Direction — cholesterol → calcium, not the reverse
_Scripts: `mr-tc-calcium.qmd`, `mr-ldlc-calcium.qmd`, `mr-calcium-tc.qmd`, `mr-calcium-ldlc.qmd` (UKB exposures → MGI-BioVU calcium)._

- Total cholesterol → calcium: β = 0.065, p = 5.8e-4 (280 SNPs)
- LDL-C → calcium: β = 0.053, p = 0.006 (232 SNPs)
- Calcium → total cholesterol: β = 0.032, p = 0.39 (null)
- Calcium → LDL-C: β = −0.001, p = 0.98 (null)

**Read:** cholesterol drives calcium; calcium does not drive cholesterol.

### 2. Drug-target specificity — it's HMGCR, not LDL-C per se
_Scripts: `drug_target_mr_analysis.qmd`, `drug_target_mr_ukb.qmd` (outcome: UKB calcium, Barton 2021)._

- HMGCR cis → calcium: β ≈ 0.185, p < 1e-20
- PCSK9 cis → calcium: β ≈ 0.004, p = 0.87 — **informatively null** (F > 2000)
- NPC1L1: imprecise / too few cis instruments to instrument reliably

**Read:** the calcium effect tracks the mevalonate pathway, not LDL-C lowering in
general. The HMGCR effect is ~3–5× the all-cholesterol effect, i.e. HMGCR-specific.

### 3. Bone & related outcomes — cholesterol lowers BMD but not fractures; vitamin D goes the wrong way
_Scripts: `drug_target_mr_bmd.qmd`, `mr-downstream-analyses.qmd`._

- HMGCR cis → Heel BMD: β ≈ −0.115 (LDL) / −0.106 (TC), p < 1e-12
- All-cholesterol → Heel BMD: β ≈ −0.051, p = 1.4e-4
- Femoral neck BMD: null (β ≈ −0.008, p = 0.90 — but underpowered, GEFOS n≈33k)
- Fractures (Dönertaş, n≈484k): **null** (β ≈ 0.0002, p = 0.73)
- Total cholesterol → 25-OH vitamin D (Revez UKB): β = **−0.137, p = 1e-29** (strong **inverse**)

**Read:** cholesterol lowers heel BMD (concordant with statin→BMD literature, and
in the *same* direction), but this does **not** translate to fracture risk, and
vitamin D moves *inversely* — so a "cholesterol → ↑vitamin D → ↑calcium" route is
refuted (wrong sign). These set up the mediation tests below.

### 4. Mediation MVMR — BMD does not carry the effect
_Script: `mvmr_analyses.qmd` (spec: `MVMR_Analysis_Specification.md`)._

- BMD → calcium (the decisive leg): β ≈ −0.005, p = 0.76 (**null**)
- HMGCR → calcium mediated by BMD (two-step): **≈ 0.5 %**
- All-LDL → calcium: total 0.039 → direct | BMD 0.027 (difference method, both non-sig)

**Read:** the chain breaks at the **second** link (BMD → calcium), not the first.
Bone loss and calcium rise are parallel HMGCR effects, consistent with serum
calcium being a homeostatically defended set-point.

### 5. Mediator second-leg screen — no bone-density site mediates; PTH untestable
_Script: `mediator_screen.qmd`._

- Heel eBMD → calcium: β ≈ −0.009, p = 0.45 (well-powered: F ≈ 197, 301 SNPs)
- Femoral neck BMD → calcium: β ≈ −0.001, p = 0.97
- Total body BMD → calcium: β ≈ −0.012, p = 0.42
- PTH (INTERVAL `prot-a-2431`, KORA `prot-c-3726_62_4`): **no cis instrument** → could not test

**Read:** no bone-density measure (heel QUS, femoral-neck or total-body DXA)
causally moves serum calcium — this closes the "maybe hip/femoral-neck mediates"
question. PTH is *untested*, not excluded (no strong cis-pQTL exists in these
datasets).

### 6. Outcome robustness — not an albumin artefact; phosphate flat → calcium-specific
_Script: `calcium_artefact_checks.qmd` (albumin & phosphate from the UKB Barton 2021 biomarker batch)._

- Albumin → calcium (positive control): β = 0.526, p = 1e-177 — assay artefact channel is open & the test is sensitive
- HMGCR → albumin: β = −0.042, p = 0.024 (**negative**)
- HMGCR → calcium, albumin-adjusted (two-step direct): ≈ 0.207 (**112 % of total**) → **not an albumin artefact**
- HMGCR → phosphate: β ≈ 0.027, p = 0.15 (**flat**)
- Phosphate → calcium: β ≈ −0.051, p = 0.37 (null)

**Read:** the calcium signal is real (albumin, if anything, opposes it). Flat
phosphate rules out PTH (would drop PO₄), vitamin-D/absorption and bone resorption
(would raise PO₄). Combined with the bone-density null, everything converges on a
**calcium-specific, plausibly renal (tubular handling)** mechanism.

---

## Key numbers at a glance

| Relationship | β (IVW/two-step) | p | Verdict |
|---|---|---|---|
| Total cholesterol → calcium | 0.065 | 5.8e-4 | causal (forward) |
| LDL-C → calcium | 0.053 | 0.006 | causal (forward) |
| Calcium → cholesterol (either) | ~0 | >0.3 | null (no reverse) |
| **HMGCR → calcium** | **0.185** | <1e-20 | strong, HMGCR-specific |
| PCSK9 → calcium | 0.004 | 0.87 | informative null |
| HMGCR → heel BMD | −0.115 | <1e-12 | lowers BMD |
| Cholesterol → fractures | ~0 | 0.73 | null |
| Cholesterol → vitamin D | −0.137 | 1e-29 | inverse (refutes +ve route) |
| **BMD → calcium** | **−0.005** | 0.76 | **null (no mediation)** |
| Femoral neck / total-body BMD → calcium | ~0 | >0.4 | null (all sites) |
| Albumin → calcium (control) | 0.526 | 1e-177 | artefact channel open |
| HMGCR → calcium \| albumin | 0.207 | — | 112 % → not artefact |
| HMGCR → phosphate | 0.027 | 0.15 | flat → calcium-specific |
| PTH → calcium | — | — | untestable (no cis instrument) |

---

## Methodological lessons (carry forward)

- **Trust the two-step (product-of-coefficients) over difference-method MVMR for
  cis exposures.** Pooling a few cis SNPs with hundreds of genome-wide mediator
  instruments contaminates the cis coefficient (bit PCSK9 in `mvmr_analyses.qmd`
  and the albumin MVMR in `calcium_artefact_checks.qmd`). The two-step never has
  to estimate a small direct effect as the difference of two large numbers.
- **Serum calcium is a homeostatically defended set-point**, not a reservoir
  readout — so graded genetic variation in *any* bone measure is expected to be
  buffered. This is the physiological reason BMD → calcium (and likely turnover →
  calcium) comes back null.
- **One-sample UKB overlap** (LDL / BMD / calcium all UKB) biases estimates toward
  the null — so positive findings here are conservative.
- **For pQTL instruments, cis strength beats discovery N.** The PTH pQTLs failed
  not because n≈1–3k is small, but because PTH has *no* strong cis-pQTL (pulsatile,
  physiologically regulated). Use a cis-pQTL (clean, avoids calcium-feedback loci)
  and colocalization instead of Egger/PRESSO when only 1 SNP is available.
- **Batch OpenGWAS `associations()` queries** (~100 variants) — a long variant list
  is silently truncated, which once collapsed an MVMR union from 528 → 80 SNPs.

---

## Open threads / next steps

| Thread | Status | Notes |
|---|---|---|
| **Renal test — urinary calcium / fractional excretion** | **not started** (highest priority) | Directly probes the leading (tubular) hypothesis. Candidate GWAS: Sun et al. 2020, 24-h urinary Ca/Mg/urate (PMID 31993563). Test HMGCR → urinary calcium. |
| **PTH via a larger pQTL** | not started | deCODE (Ferkingstad 2021, SomaScan, n≈35k, non-UKB) preferred; UKB-PPP Olink is bigger but overlaps UKB calcium. Needed to *test* (not just assume) the PTH route. |
| **Bone-turnover markers (CTX, P1NP, osteocalcin)** | blocked on data | Mechanistically motivated (net-mineral-balance model) but no public GWAS wired in; add a local loader if summary stats obtained. |
| **HDL-C / triglyceride associations** | not started | From `TODO.md`. |
| **PCSK9 pQTL MR (UKB-PPP)** | not started | From `TODO.md` — strengthen the PCSK9-null leg with protein-level instruments. |
| **Colocalization of HMGCR signals (LDL-C / BMD / eQTL)** | not started | From `TODO.md` — H4 vs H3 to defend cis-instrument validity. |
| **HMGCR eQTL-based MVMR across tissues** | stub only | `mvmr_analyses.qmd` Part 4 is a guarded stub pending eQTL summary stats. |
| **Sample-overlap correction method** | not decided | From `TODO.md`. |

---

## Script index

| Script | Produces / role |
|---|---|
| `mr-{tc,ldlc}-calcium.qmd`, `mr-calcium-{tc,ldlc}.qmd` | Bidirectional cholesterol ↔ calcium MR (direction) |
| `drug_target_mr_analysis.qmd`, `drug_target_mr_ukb.qmd` | Drug-target (HMGCR/PCSK9/NPC1L1) → calcium |
| `drug_target_mr_bmd.qmd` | Drug-target → BMD / fractures / vitamin D |
| `mr-downstream-analyses.qmd` | Cholesterol → vitamin D, BMD, fractures (secondary outcomes) |
| `mvmr_analyses.qmd` | Mediation MVMR: does BMD carry cholesterol → calcium? (spec: `MVMR_Analysis_Specification.md`) |
| `mediator_screen.qmd` | Second-leg screen: which mediators affect calcium? |
| `calcium_artefact_checks.qmd` | Outcome robustness: albumin (artefact) & phosphate (co-regulation) |
| `summary-tables.qmd` | Collates outputs for figures/tables |

See `README.md` for execution order, datasets, and software/reproducibility notes.

---

## Changelog

- **2026-07-11** — Initial version. Captures the arc through the mediation MVMR,
  mediator second-leg screen, and albumin/phosphate outcome-robustness checks:
  bone-density mediation and albumin artefact excluded; phosphate flat; mechanism
  narrowed to calcium-specific / renal. PTH untested; renal test is the next step.
