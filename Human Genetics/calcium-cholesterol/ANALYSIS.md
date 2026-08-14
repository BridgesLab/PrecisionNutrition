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

### 3b. Robustness of cholesterol → BMD to sample overlap and correlated pleiotropy
_Scripts: `robust_mr_prep.qmd`, `robust_mr_apss.qmd`, `robust_mr_cause.qmd`,
`robust_mr_summary.qmd`. Rationale: `ROBUST_MR.md`._

The heel-BMD result above rests on a UK Biobank exposure against a UK Biobank
outcome. Two genome-wide methods were applied that model this explicitly:
**MR-APSS** (sample structure via bivariate LDSC intercepts, plus a winner's-curse
correction licensing a relaxed 5e-5 instrument threshold) and **CAUSE**
(shared-factor model for correlated horizontal pleiotropy). Each was run on an
overlapping exposure (GLGC 2021, N≈1.32M) and an overlap-free one (GLGC 2013,
entirely pre-UK-Biobank), so the overlap correction is *measured* rather than
assumed.

**Overlap is not reliably measurable, and does not matter.** Three independent
estimators, from three model families:

| | MR-APSS C₁₂ | CAUSE ρ | MRBEE R_xy |
|---|---|---|---|
| GLGC 2021 (overlapping) | +0.0226 (SE 0.0122) | −0.0109 | −0.0045 |
| GLGC 2013 (pre-UKB) | −0.0074 (SE 0.0095) | −0.0146 | −0.0023 |

They disagree on sign, and CAUSE's is *larger* in the arm with no overlap by
construction — so none of them is resolving overlap precisely. What they agree on
is magnitude: all are within ±0.023 of zero. Switching MR-APSS's correction on
moves the estimate **0.17 SE** (0.07 SE in the clean arm). Overlap is immaterial
here, and we should not claim to have detected it.

**Threshold matters more than method class.** Comparing like-for-like at a fixed
instrument threshold, the two assumption classes agree closely:

| Method | @ 5e-8 | @ 5e-5 |
|---|---|---|
| IVW | −0.067 | — |
| MR-RAPS (l2) | −0.084 | −0.056 |
| MR-RAPS (tukey, robust) | −0.062 | −0.043 |
| MRBEE | −0.072 | −0.051 |
| MR-APSS | — | −0.056 |
| CAUSE γ | — | −0.040 |

(overlapping arm; the overlap-free arm shows the same pattern at ~60 % the
magnitude). MR-RAPS(l2) and MR-APSS land on −0.056 and −0.056 at the same
threshold despite opposite assumptions about correlated pleiotropy. **The apparent
"more conservative method, smaller estimate" gradient is largely a threshold
effect**, not evidence about pleiotropy.

**Two corrections that nearly cancel.** Weak-instrument correction pushes away
from the null (IVW −0.067 → RAPS-l2 −0.084 on identical instruments); robust loss
pushes back (−0.084 → −0.062). IVW lands between them by coincidence, not validity.

**Pleiotropy is heavy but diffuse.** MR-RAPS estimates τ̂² ≈ 2×10⁻⁵, i.e. a
pleiotropy SD of ~0.005 — roughly 3× the outcome standard error. Yet MRBEE's
iterative outlier test flags only 1–3 % of instruments (7/443, 12/948). Those are
consistent only if pleiotropy is spread thinly across many variants rather than
concentrated in outliers, which favours a random-effects treatment (RAPS) over
deletion (MRBEE) — and explains why MRBEE tracks RAPS-l2 rather than RAPS-tukey.

**IVW's precision was indefensible.** Its SE of 0.0056 is contradicted by RAPS,
MRBEE and MR-APSS independently, all of which put it at 0.017–0.027. Three model
families agree the conventional interval was ~4× too narrow.

**CAUSE cannot discriminate.** Sharing vs causal: Δelpd = −1.54 (SE 1.53, p = 0.16)
and −0.40 (SE 1.28, p = 0.38). Null vs causal is *also* non-significant, so no
model wins — every ELPD difference is smaller than its own standard error.
Pseudo-BMA weights put 0.74 of the predictive support on the causal model in the
overlapping arm but only 0.40 in the underpowered one. Read Bayesianly and
conditional on the causal model, P(γ < 0) ≈ 0.98 and ≈ 0.92.

**Read:** cholesterol → heel BMD is negative and non-null across five estimators
spanning both assumption classes, at every threshold, in both exposure arms. The
effect is real. What changes is precision, not direction: the honest interval is
roughly −0.09 to −0.01 rather than IVW's −0.078 to −0.056.

Neither of the two things this arm was built to test turns out to be the
explanation for anything. Sample overlap is immaterial (0.17 SE) and not reliably
detectable. Correlated pleiotropy cannot be *excluded* — not because CAUSE finds
evidence for it (q ≈ 0.02–0.03) but because CAUSE lacks the power to exclude
anything. The variation across methods is driven by instrument threshold and by
how each handles diffuse balanced pleiotropy.

This is why the HMGCR cis result (§2, §7) carries the argumentative weight: it is
a design these genome-wide polygenic methods cannot address, and cannot undermine.

**Attenuation in the clean arm is weak instruments, not overlap.** GLGC 2013 has
75 instruments at F = 61 versus GLGC 2021's 443 at F = 129, and weak instruments
bias two-sample MR toward the null — a simpler explanation than overlap inflation,
and the one both the C = I comparison and the RAPS threshold contrast support.

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

### 7. HMGCR expression → BMD (blood eQTL) — one solid result; total-expression MR otherwise not viable
_Script: `hmgcr_eqtl_mvmr.qmd`. Thread stopped here by design (see limitation)._

- **HMGCR expression (whole blood, eQTLGen) → heel BMD: β = −0.044, se = 0.011,
  p = 4.4e-5** (single cis-eQTL rs6453133, F = 224, Wald ratio). Direction is as
  predicted (↑HMGCR expression → ↓BMD, matching the cholesterol cis-MR of −0.115).
- Because whole blood proxies the **myeloid/osteoclast lineage**, this supports the
  **cell-autonomous osteoclast model** (mevalonate → GGPP prenylation in the
  osteoclast precursor) — a solid, correctly-signed, well-instrumented result.

**Limitation — why total-expression eQTL MR is not the tool for HMGCR (thread closed):**
- The **functional HMGCR variants are null for total mRNA** in GTEx: rs3846662
  (exon-13 splice variant) p = 0.78 (blood) / 0.93 (liver); rs12916 p ≈ 0.46;
  rs6453133 p ≈ 0.65. **HMGCR is regulated by *splicing*, not expression level**, so
  a gene-expression (`ge`) eQTL is blind to the causal mechanism.
- GTEx per-tissue power is far too low regardless (best HMGCR `ge` p ≈ 1e-3;
  kidney n = 73), so the multi-tissue calcium/kidney/liver screen and colocalization
  cannot run.
- The eQTLGen blood signal exists only because n ≈ 31k detects a tiny total-
  expression effect; its lead (rs6453133) is multiallelic and was dropped by the
  calcium GWAS's biallelic QC, so **blood → calcium could not be tested**.

**Future avenues (if ever revisited):** (a) **HMGCR splice-QTL** (eQTL Catalogue
`leafcutter`/`txrev`; the rs3846662 exon-13 event) — the functionally correct,
strong, biallelic exposure; (b) eQTLGen-full blood `ge` for a well-powered (if
biologically small) total-expression instrument that dodges the multiallelic SNP.
Neither is critical: the **drug-target cis-MR already provides a valid HMGCR
instrument** (rs3846662/rs12916), so tissue-of-action is the only thing the eQTL
route would add — and that is limited by cis-eQTL sharing across tissues anyway.

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
| **HMGCR expression (blood eQTL) → BMD** | **−0.044** | 4.4e-5 | solid; ↓BMD, osteoclast-lineage |
| HMGCR expression (blood eQTL) → calcium | — | — | untestable (multiallelic lead absent from calcium GWAS) |
| LDL-C → heel BMD, MR-APSS (GLGC 2021) | −0.056 | 0.034 | survives overlap + pleiotropy modelling |
| LDL-C → heel BMD, MR-APSS (GLGC 2013, no overlap) | −0.037 | 0.045 | replicates in an overlap-free arm |
| LDL-C → heel BMD, CAUSE γ (GLGC 2021) | −0.040 | 0.16 | directionally consistent; model comparison inconclusive |
| LDL-C → heel BMD, MR-RAPS l2 @5e-5 (GLGC 2021) | −0.056 | 0.007 | matches MR-APSS exactly, opposite assumptions |
| LDL-C → heel BMD, MR-RAPS tukey @5e-5 | −0.043 | 0.010 | robust loss; outliers inflate the l2 fit |
| LDL-C → heel BMD, MRBEE @5e-5 (GLGC 2021) | −0.051 | 0.009 | overlap + weak-IV corrected |
| Sample-overlap effect on the estimate (MR-APSS C vs C=I) | 0.004 | — | **0.17 SE — immaterial** |
| Horizontal pleiotropy SD (MR-RAPS τ̂) | 0.005 | — | ~3× the outcome SE — heavy but balanced |

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
- **One-sample UKB overlap biases toward the null — but the magnitude here is
  negligible, and we have now measured it rather than assumed it** (§3b). MR-APSS
  puts the cross-trait LDSC intercept at C₁₂ = 0.023 for GLGC 2021 × UKB eBMD
  versus ~0 for a pre-UKB exposure, confirming shared samples; yet switching the
  correction on moves the estimate 0.17 SE, and CAUSE's ρ is unchanged between
  arms. Direction of the old assumption was right, magnitude was not worth the
  worry. Do not spend further effort on overlap for this exposure/outcome pair.
- **Attenuation between a big overlapping GWAS and a small clean one is usually
  weak instruments, not overlap.** GLGC 2013 gives a smaller effect than GLGC 2021,
  which looks like overlap inflation until you notice F = 61 vs 129 and that the
  explicit overlap correction does nothing.
- **Compare methods at a fixed instrument threshold, or you are measuring the
  threshold.** Across five estimators the p-value threshold moved the estimate
  more than the choice of method did: everything sits near −0.07 at 5e-8 and near
  −0.05 at 5e-5, and MR-RAPS(l2) and MR-APSS agree to three decimals at 5e-5
  despite opposite assumptions about correlated pleiotropy. A "more conservative
  method gives a smaller estimate" narrative is easy to construct by accident.
- **Distinguish diffuse from concentrated pleiotropy before choosing a method.**
  Here τ̂ ≈ 0.005 (≈3× the outcome SE) while MRBEE's outlier test flags only 1–3 %
  of instruments — pleiotropy spread thinly across many variants. Random-effects
  methods (RAPS) suit that; outlier-deletion methods (MRBEE, MR-PRESSO) are built
  for the opposite case and will under-correct.
- **Genome-wide polygenic MR methods cannot adjudicate a cis design.** CAUSE and
  MR-APSS need thousands of instruments and a polygenic background model; they say
  nothing about the HMGCR cis result, and a null from them would not undermine it.
  Keep the two arms rhetorically separate.
- **Verify GWAS-Catalog column semantics, never trust the header.** In GCST006979
  (Morris 2019 eBMD) `variant_id` is `chr:pos:ref:alt` and the rsIDs live in a
  column called `snp.1`; the column labelled `n` is BOLT-LMM's second p-value, not
  a sample size. Mapping the obvious names silently produced a zero-row dataset.
  `read_gwas()` now cross-checks any p-value column against the one implied by
  beta/se, which catches this class of error immediately.
- **For pQTL instruments, cis strength beats discovery N.** The PTH pQTLs failed
  not because n≈1–3k is small, but because PTH has *no* strong cis-pQTL (pulsatile,
  physiologically regulated). Use a cis-pQTL (clean, avoids calcium-feedback loci)
  and colocalization instead of Egger/PRESSO when only 1 SNP is available.
- **Batch OpenGWAS `associations()` queries** (~100 variants) — a long variant list
  is silently truncated, which once collapsed an MVMR union from 528 → 80 SNPs.
- **Check the molecular phenotype before an expression-eQTL MR.** HMGCR is
  regulated by *splicing*, so its functional variants (rs3846662, rs12916) are null
  for total mRNA (`ge`) — a gene-expression eQTL MR is blind to it. Use the QTL
  type that matches the biology (sQTL/txrev here), and don't rely on GTEx for a
  gene whose per-tissue eQTL is weak (small n; best HMGCR `ge` p ≈ 1e-3).

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
| **HMGCR eQTL-based MR across tissues** | **CLOSED — limitation recorded** (`hmgcr_eqtl_mvmr.qmd`, evidence-chain §7) | Kept the one solid result (blood eQTL → BMD). Total-expression eQTL MR not viable for HMGCR (splice-regulated; functional variants null for `ge`; GTEx underpowered). Reopen only via splice-QTL (`leafcutter`/rs3846662) or eQTLGen-full — not critical since the drug-target cis-MR already gives a valid HMGCR instrument. |
| **Sample-overlap correction method** | **CLOSED — resolved 2026-08-11** (`robust_mr_*.qmd`, evidence-chain §3b) | Settled on MR-APSS (LDSC-based C matrix + winner's-curse correction) and CAUSE (shared-factor model, robust to overlap via ρ), each run on an overlapping and an overlap-free exposure arm. Overlap is detectable (C₁₂ = 0.023) but immaterial (0.17 SE). No further work needed for this exposure/outcome pair. |
| **MRAID** | deferred by design | Excluded from §3b because it is a strict two-sample method and UKB sits on both sides. Appropriate under a genuinely non-overlapping design — **UKB exposure → MGI/BioVU outcome**, the same architecture already used for calcium in `mr-tc-calcium.qmd`. Blocked on an MGI or BioVU BMD/fracture phenotype. |
| **Winner's-curse correction, quantified** | not started | §3b ran MR-APSS with `Cor.SelectionBias = TRUE` throughout, so we know what the *overlap* correction does (nothing) but not what the *selection* correction does. A third fit with `Cor.SelectionBias = FALSE` would separate the two, ~2 min of compute. |

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
| `robust_mr_prep.qmd` | Genome-wide sumstats QC/harmonisation for the robust MR arm + IVW sanity anchor |
| `robust_mr_apss.qmd` | MR-APSS: sample structure (C), correlated pleiotropy (Ω), winner's curse |
| `robust_mr_cause.qmd` | CAUSE: shared-factor model + Bayesian reading of γ |
| `robust_mr_raps.qmd` | MR-RAPS: random-effects pleiotropy (τ²) + weak-instrument profile score |
| `robust_mr_mrbee.qmd` | MRBEE: bias-corrected estimating equation; overlap via R_xy |
| `robust_mr_summary.qmd` | Reconciles robust vs conventional estimates (rationale: `ROBUST_MR.md`) |
| `summary-tables.qmd` | Collates outputs for figures/tables |

See `README.md` for execution order, datasets, and software/reproducibility notes.

---

## Changelog

- **2026-07-11** — Initial version. Captures the arc through the mediation MVMR,
  mediator second-leg screen, and albumin/phosphate outcome-robustness checks:
  bone-density mediation and albumin artefact excluded; phosphate flat; mechanism
  narrowed to calcium-specific / renal. PTH untested; renal test is the next step.
- **2026-07-11** — Scaffolded `hmgcr_eqtl_mvmr.qmd`: multi-tissue HMGCR
  expression → calcium & BMD eQTL MR + colocalization, including a monocyte /
  whole-blood osteoclast-lineage proxy (bone/marrow eQTL is absent from GTEx).
  Runs for eQTLGen whole blood; awaits GTEx/BLUEPRINT accessions or local files.
- **2026-07-13** — **Closed the eQTL thread** (evidence-chain §7). Kept the solid
  result: HMGCR expression (blood) → BMD, β = −0.044, p = 4.4e-5, correct sign,
  osteoclast-lineage. Recorded the limitation: HMGCR is splice-regulated, so
  total-expression (`ge`) eQTL MR is not viable (functional variants null for `ge`;
  GTEx underpowered; eQTLGen lead multiallelic & absent from the calcium GWAS).
  Future avenues noted (splice-QTL / eQTLGen-full); not critical.
- **2026-08-11** — **Closed the sample-overlap thread** (evidence-chain §3b). Added
  an overlap- and pleiotropy-robust arm for cholesterol → heel BMD: MR-APSS and
  CAUSE, each on an overlapping (GLGC 2021) and an overlap-free (GLGC 2013,
  pre-UKB) exposure. Overlap is detectable (C₁₂ = 0.023 vs ~0) but immaterial
  (0.17 SE). The effect survives at −0.056 (MR-APSS, p = 0.034) and −0.037 in the
  clean arm; CAUSE is directionally consistent (γ = −0.040) but its model
  comparison is inconclusive in both arms — it cannot separate causal from sharing
  *or* from null. MRAID excluded by design (two-sample assumption violated);
  deferred to a UKB → MGI/BioVU design. Rationale and full method matrix in
  `ROBUST_MR.md`.
- **2026-08-12** — Extended §3b with the InSIDE-assumption class: **MR-RAPS**
  (`robust_mr_raps.qmd`) and **MRBEE** (`robust_mr_mrbee.qmd`), both run at 5e-8
  and 5e-5 on the same harmonised substrate. Two revisions to the previous
  entry: (1) the "overlap is detectable" claim is **withdrawn** — a third
  estimator (MRBEE R_xy = −0.0045) disagrees in sign with MR-APSS C₁₂ (+0.023),
  so overlap is immaterial but not reliably measurable; (2) the apparent
  method-conservatism gradient is largely a **threshold** effect — at a fixed 5e-5
  threshold MR-RAPS(l2) and MR-APSS both give −0.056 despite opposite assumptions
  about correlated pleiotropy. New finding: pleiotropy is heavy (τ̂ ≈ 0.005, ~3×
  the outcome SE) but diffuse, so IVW's SE was ~4× too narrow and outlier-deletion
  methods under-correct.
