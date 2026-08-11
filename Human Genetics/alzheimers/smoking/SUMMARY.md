# SUMMARY — Decomposing the Smoking → Alzheimer's MR Signal

*Pipeline run 2026-06-04. Exposure: cigarettes/day (GSCAN/Liu 2019, `ieu-b-142`). Primary
outcome: AD (Bellenguez 2022, `ebi-a-GCST90027158`). **Canonical baseline = the 20-SNP set
with the CLU/MINDY2 horizontal-pleiotropy outliers removed (Layer 0b)**; Layers 1–5 build on
it. All estimates on the log-odds scale per SD of cigarettes/day unless noted. Full decision
trail in [`results/decisions_log.md`](results/decisions_log.md).*

## Answers to the four review questions

1. **Canonical baseline is now the outlier-removed estimate** — IVW b = −0.113, SE = 0.032,
   **p = 3.4×10⁻⁴** (20 SNPs), and Layers 1–5 read the pruned set / pruned instrument file.
2. **Any other clusters once CLU/MINDY2 are gone?** No. MR-Clust on the full set still
   isolates CLU/MINDY2 as a distinct *risk* cluster (θ = +1.24); on the pruned set the
   remaining 20 SNPs form a **single protective cluster with no further substructure**.
3. **Did the relaxed criteria add SNPs / new targets?** Yes for instruments — 15q25
   CHRNA5/A3/B4 went 0→31 SNPs and 8p11 CHRNA6/B3 1→6 (pooled nAChR 3→43). No *new genes*:
   Layer 6 is cis-restricted to the receptor panel, so it newly *instruments* known loci
   rather than nominating off-panel targets (a genome-wide relaxed + druggable-genome scan
   would be needed for that — not yet run).
4. **Is 15q25 one LD block — can we separate CHRNA5/A3/B4?** The block carries **7
   independent signals** (r²<0.01) but they span the three overlapping subunit genes
   (max pairwise |r| = 0.56), so MR gives only a **locus-level** effect — subunits are
   **not separable** without subunit-specific brain QTL + conditional/coloc. And CHRNA4 is
   **not** significantly stronger than 15q25 (Δb = −0.15, z = −0.97, **p = 0.33**): its
   larger point estimate just has a wide CI (4 SNPs); 15q25's smaller p reflects more
   instruments, not a bigger effect.

## Headline

When smoking instruments are selected **at the nAChR receptor loci using drug-target-MR
conventions** (relaxed p<1e-5, correlated cis-SNPs r²<0.3, LD-aware IVW; Layer 6), the
protective smoking→AD effect is **strong and significant at the receptor loci**: pooled
across the nAChR panel, IVW b = −0.129 (OR 0.88, p = 9×10⁻⁴, 43 cis-SNPs), driven by the
**15q25 CHRNA5/A3/B4** cluster and **CHRNA4**. This is a large gain over the strict
genome-wide instrument set (3 SNPs, p = 0.06) and over the diluted composite signal — which
itself becomes significant once two AD-direct horizontal-pleiotropy outliers (**CLU**,
**MINDY2**) are removed (IVW b = −0.113, p = 3×10⁻⁴).

Critically, with **locus-based binning** (assign SNPs by panel-gene cis-window, not literal
nearest gene), the mechanism-stratified MR now **localises** the protection: the nAChR bin is
significantly protective (b = −0.149, p = 4×10⁻⁴, 4 SNPs) while the "other pleiotropic" bin
collapses to **null** (b ≈ 0, p = 0.95). The protection is therefore **mechanism-specific
(nicotinic), not diffuse genome-wide pleiotropy**.

What is **still not demonstrated** is the *behavior-independent, druggable* molecular effect:
the only receptor gene with a (blood) cis-eQTL, **CHRNB2**, shows a null expression→AD
estimate (p = 0.79), and the brain-expressed receptor genes have no blood-eQTL instrument.
So the evidence now strongly localises the protection to the **nAChR receptor loci**, but
isolating receptor *function* (the truly druggable quantity) still requires **brain
eQTL/pQTL + colocalization**. No drug target is formally nominated yet, but CHRNA5/A3/B4 and
CHRNA4 are now the clear priority targets to pursue with brain QTL.

⚠️ **Read the headline magnitude with a major caveat.** The **Selection/collider-bias section**
(Layers 4b–4e) now finds the protective estimate is **most plausibly a survival/competing-risk
collider artefact**: a properly identified index-event correction (2,346 LD-clumped SNPs,
b_SH = −0.944) reduces the pooled effect from **−0.129 to +0.005 (p=0.91)**, and the collider
structure quantitatively predicts the observed magnitude (−0.134 predicted vs −0.129 observed).
Proxy-ascertainment bias is *not* the mechanism (the effect persists in clinical-only AD) —
survival selection is. **The protective effect below should not be treated as a causal,
druggable signal** pending an incident / younger-onset AD replication.

## Layer-by-layer

### L0 — Composite signal (the baseline to decompose)
Direction is protective and robust across estimators; magnitude modest; substantial
heterogeneity (I² ≈ 75%), and the more outlier-robust estimators are *more* protective —
consistent with a few risk-direction pleiotropic SNPs diluting the IVW mean.

| Method | b (log OR) | OR | p |
|---|---|---|---|
| IVW (multiplicative RE) | −0.065 | 0.94 | 0.296 |
| IVW (fixed effects) | −0.065 | 0.94 | 0.037 |
| MR-RAPS | −0.098 | 0.91 | 0.050 |
| MR-Egger | −0.208 | 0.81 | 0.059 |
| Weighted median | −0.134 | 0.87 | **0.0017** |
| Weighted mode | −0.130 | 0.88 | 0.647 |
| MR-PRESSO (outlier-corrected) | −0.112 | 0.89 | **0.0015** |

Egger intercept p = 0.108 (no strong *directional* pleiotropy). → **Proceed to decompose.**

**Pleiotropy-outlier sensitivity (L0b).** Two SNPs — rs73229090 (**CLU**, a canonical AD
gene) and rs632811 (**MINDY2**) — behave as uncorrelated horizontal pleiotropy: MR-Clust
isolates them into a separate *risk* cluster (L1), they dominate leave-one-out, and they are
the MR-PRESSO outliers. Their effect is direct-to-AD, not via smoking. Removing them moves
IVW from b = −0.065 (p = 0.30) to **b = −0.113, p = 3.4×10⁻⁴** (weighted median b = −0.134,
p = 0.001) — i.e. they were *masking* a significant protective effect.

| Set | IVW b | p | Weighted-median b | p |
|---|---|---|---|---|
| all 22 SNPs | −0.065 | 0.296 | −0.134 | 0.0017 |
| **drop CLU + MINDY2 (20 SNPs)** | **−0.113** | **0.00034** | −0.134 | 0.0011 |

### L1 — Discovery (MR-Clust), run on two sets
**Full 22-SNP set (discovery):** one dominant protective cluster (20 SNPs, θ ≈ −0.07) plus a
distinct **risk** cluster of exactly 2 SNPs — rs73229090 (**CLU**) and rs632811 (**MINDY2**),
θ ≈ +1.24. This *is* the formal identification of the two pleiotropy outliers.

**Pruned canonical set (CLU/MINDY2 removed):** the remaining 20 SNPs collapse into a
**single protective cluster with no further substructure** — there are **no other relevant
clusters**. Note MR-Clust works in *effect-ratio* space, so it does not separate the nAChR
SNPs as their own cluster; receptor-locus specificity instead comes from the **locus-based
mechanism stratification (L2)** and the relaxed cis analysis (L6).

### L2 — Mechanism-stratified MR (pruned set, **locus-based binning**)
SNPs are binned by **cis-window membership of the panel genes** (with a nearest-gene fallback),
not literal nearest gene — so the 15q25 cluster lead (nearest *HYKK*) is correctly assigned to
nAChR, and the chr19 locus (nearest a pseudogene/*EGLN2*) to nicotine_metabolism. This is the
decisive decomposition:

| Bin | nSNP | IVW b | IVW p | Weighted-median b | WM p |
|---|---|---|---|---|---|
| **nAChR pharmacodynamic** (15q25, CHRNA4, CHRNB2, CHRNB3) | 4 | **−0.149** | **0.0004** | **−0.145** | **0.001** |
| nicotine_metabolism (CYP2A6 region) | 2 | −0.117 | 0.122 | – | – |
| reward/behavioral (DRD2, DBH) | 2 | −0.131 | 0.378 | – | – |
| other pleiotropic | 12 | **−0.005** | **0.946** | 0.048 | 0.625 |

**The protection localises.** Once the SNPs are correctly assigned, the **nAChR bin is
significantly protective (p = 4×10⁻⁴)** while the **"other" bin collapses to null (b ≈ 0,
p = 0.95)**. The earlier appearance of a protective "other" bin was an artefact of mis-binned
receptor/metabolism SNPs (the strongly protective 15q25 and chr19 loci) sitting in it. The
nicotine-metabolism bin trends protective too. This is the cleanest evidence that the signal
is **mechanism-specific (nicotinic), not diffuse genome-wide pleiotropy**.

### L3 — Druggable direct effect: blocked by blood-eQTL coverage
Of the 8-gene receptor panel, only **CHRNB2** (`eqtl-a-ENSG00000160716`) exists as a blood
eQTL dataset in OpenGWAS; the other seven are brain-predominant and have no blood cis-eQTL
instrument there. For CHRNB2 the behavior-independent test is **null**:

- cis-MR (expression → AD): b = 0.002, p = 0.96 (1 cis SNP)
- MVMR (CHRNB2 expression | smoking): receptor b ≈ 0 (p = 0.998); smoking retains b = −0.12 (p = 0.11)
- colocalization with AD: **PP.H4 = 0.004** (no shared causal variant; PP.H3 vs H4 favors distinct signals)

→ **0 genes nominated** (threshold: coloc PP.H4 ≥ 0.8 *and* cis-MR p < 0.05). See
[`results/targets.tsv`](results/targets.tsv).

### L4 — Robustness
- **Survival/selection collider:** the CHRNB2 instrument is ~null on parental lifespan
  (p = 0.43), lung cancer (p = 0.77) and CAD (p = 0.059) — i.e. the behavioral/mortality
  axis looked amputated **for the one gene testable**. ⚠️ **Superseded by L4b:** CHRNB2 is the
  ≈null *contributor*; when the **effect-driving** 15q25/nAChR instruments are tested they are
  emphatically **not** selection-clean (15q25 → lung cancer p=10⁻⁵¹). Do not read this CHRNB2
  result as evidence the signal is collider-free — see the Selection/collider-bias section.
- **Measured pleiotropy (MVMR on BMI/LDL/EA):** the smoking→AD effect *survives and
  strengthens* after conditioning — direct smoking b = −0.092 (p = 0.029); it is not
  explained by BMI, LDL, or educational attainment.
- **Correlated/directional pleiotropy:** MR-Egger-intercept proxy p = 0.108 (no strong
  signal; full CAUSE needs genome-wide stats — see README).
- **Ascertainment (Bellenguez direct-case vs "Wightman" proxy/GWAX):** direct-case b = −0.065;
  the proxy estimate is near-null but on a different effect scale, so this comparison was only
  weakly informative as run. ⚠️ **Superseded by L4c:** the "Wightman" ID (`ieu-b-5067`) was
  mislabelled (actually Woolf 2022, 954 cases); the proper proxy→clinical gradient (Kunkle /
  Lambert / Bellenguez / Schwartzentruber-GWAX) on a harmonised scale is in the
  Selection/collider-bias section and shows the effect **persists** in clinical-only AD.

### L6 — Relaxed cis drug-target MR at nAChR loci (LD-aware): the key gain
Selecting smoking instruments *within each receptor cis window* at drug-target thresholds
(p<1e-5, r²<0.3, LD-aware IVW/Egger) recovers far more receptor-locus signal — most
importantly the **15q25 CHRNA5/A3/B4 cluster**, which strict genome-wide clumping + nearest-
gene binning had pushed into "other".

| Locus | nSNP | IVW b | IVW p | Egger b | Egger p |
|---|---|---|---|---|---|
| **CHRNA5/A3/B4 (15q25)** | 31 | −0.105 | **0.023** | −0.133 | **0.006** |
| CHRNA4 | 4 | −0.256 | 0.087 | −0.438 | 0.129 |
| CHRNB2 | 2 | −0.326 | 0.095 | – | – |
| CHRNA6/B3 (8p11) | 6 | −0.092 | 0.554 | −0.327 | 0.213 |
| CHRNA7 (15q13) | 0 | – | – | – | – (no smoking-cis SNP; α7 ≠ heaviness locus) |
| **Pooled nAChR** | **43** | **−0.129** | **0.00087** | **−0.164** | **0.00013** |

Egger ≈ IVW at the pooled level ⇒ no major directional pleiotropy *at the receptor loci*.
Strict (3 SNPs, p = 0.063) → relaxed (43 SNPs, **p = 0.0009**).

**SNP gain (what the relaxed criteria bought).** Per-locus instrument counts, strict (Layer
2) → relaxed cis (Layer 6): 15q25 **0 → 31**, 8p11 CHRNA6/B3 **1 → 6**, CHRNA4 1 → 4, CHRNB2
1 → 2, CHRNA7 0 → 0. The big gain is the previously mis-binned 15q25 cluster.

**15q25 LD structure (can we name the subunit?).** The block holds **7 independent signals**
(strict r²<0.01 clumping) but they span the three overlapping subunit genes (max pairwise
|r| = 0.56 among the relaxed instruments). The 15q25 estimate is therefore a **locus-level**
effect — CHRNA5/A3/B4 **cannot be separated** by single-trait cis-MR; that needs subunit-
specific brain eQTL/pQTL with conditional/multi-signal colocalization. CHRNA4 (chr20) and
CHRNB2 (chr1) are on other chromosomes and so *are* distinguishable from 15q25.

**Is CHRNA4 stronger than 15q25?** No. Despite a larger point estimate (b = −0.26 vs −0.10),
the difference is not significant (Δb = −0.15, z = −0.97, **p = 0.33**); CHRNA4's CI is wide
(4 SNPs), and 15q25's smaller p reflects its 31-SNP instrument count, not a larger effect.

**Molecular (cis-eQTL) exposure, relaxed:** only **CHRNB2** has a blood eQTL (now 12 cis-SNPs
vs 1 before); expression→AD remains **null** (IVW b = −0.011, p = 0.79). The brain-expressed
receptor genes still have no blood-eQTL instrument.

## Interpretation

The picture after the relaxed drug-target layer is markedly stronger than the composite
signal suggested. The protective smoking→AD effect is **robust, not driven by measured
confounders (L0, L4), not driven by the CLU/MINDY2 pleiotropy outliers (L0b), and strongly
localised to the nAChR receptor loci** — pooled receptor-locus IVW p = 9×10⁻⁴, led by
CHRNA5/A3/B4 and CHRNA4 (L6). The one *molecularly* testable receptor instrument (CHRNB2) is
clean of the survival collider — but that is the ≈null contributor; the **effect-driving 15q25
instruments are not** (L4b). ⚠️ **This "robust" reading is now substantially undercut:** the
signal survives proxy de-selection (L4c), but an identified genome-wide index-event correction
removes it entirely (−0.129 → +0.005, L4e), both collider arms are present (L4b), and the
collider structure predicts the observed magnitude. The protective effect is **most plausibly a
survival-collider artefact** — see the Selection/collider-bias section for the full verdict
before relying on anything in this paragraph.

What remains **un-demonstrated** is the decisive *druggable* quantity: a behavior-
*independent* receptor **function/expression** → AD effect that **colocalizes** at a locus
(L3). The smoking-cis instruments sharpen *localisation and power* but still tag the
behavioral axis at the receptor locus — they do not by themselves prove a function-mediated
mechanism. The single molecular test available (CHRNB2 blood eQTL) is null, and brain
eQTL/pQTL for the priority genes were not reachable in this run.

Net: this is now a **strong, well-localised receptor-locus signal with named priority
targets (CHRNA5/A3/B4, CHRNA4)** — upgraded from the earlier "promising but diffuse" read —
but not yet a formally nominated drug target, because the function-level colocalization
evidence still needs brain QTL.

## Selection / collider-bias decomposition (Layers 4b–4e, added 2026-06-16)

The open inferential threat was **selection / survival collider bias** (not confounding or
diffuse pleiotropy, already handled by L0b/L2/L4-MVMR), made live by two facts: the primary
outcome **Bellenguez is proxy-majority** (≈39,106 clinical vs ≈46,828 UKB by-proxy cases), and
the signal-carrying locus **15q25** is the canonical lung-cancer / COPD / shortened-lifespan
locus. The prior L4 selection check tested only **CHRNB2** — the ≈null contributor, not the
effect-driving locus. Layers 4b–4e test the **effect-driving** instruments directly.
Equivalence bounds were **pre-registered** (δ=0.02 primary, δ=0.05 sensitivity; APPROACH.md
§13) *before* running. Two original config IDs were found **mislabelled** during composition
audit (`ieu-a-1001` is Okbay EA, not lifespan; `ieu-b-5067` is Woolf AD with 954 cases, not
Wightman) and replaced with verified accessions.

| Layer | Test | Result | Reading |
|---|---|---|---|
| **4b (A, X→S arm)** | 15q25 / nAChR instruments vs selection axis (TOST) | 15q25 → lung cancer b=1.36 **p=1×10⁻⁵¹**, parental lifespan b=0.13 **p=1×10⁻²¹**, COPD/FEV1-FVC **p=5×10⁻¹¹**; **0/5 sets equivalence-null** | **X→S arm present** — the effect-driving instruments are NOT selection-clean |
| **4b (Y→S arm)** | AD liability (clinical Kunkle instruments) → parental lifespan | AD → SHORTER life, all b=**+0.043 p=5×10⁻³¹**; **excl-APOE b=+0.020 p=0.003** (outcome oriented +ve=shorter, verified via APOE β=+0.057 & CHRNA5 β=+0.025) | **Y→S arm present**, incl. the non-APOE component relevant to the receptor loci (modest, b≈0.02) → **both collider arms confirmed** |
| **4c (B)** | Outcome ascertainment swap (proxy→clinical) | 15q25 clinical (Kunkle −0.136 / Lambert −0.123) vs proxy-majority Bellenguez −0.105, **p_diff=0.74**; pooled clinical −0.117 vs −0.129, p_diff=0.85 | effect **does NOT attenuate** in de-selected clinical AD → **proxy-ascertainment bias specifically is excluded** |
| **4d (E)** | Outcome direction negative controls | lung cancer→AD OR=0.93 **p=0.0036**, CAD→AD OR=0.90 **p=2×10⁻⁷** (both spuriously protective) | the outcome **carries a survival/competing-risk signature** independent of any nAChR claim |
| **4e (C)** | Exposure∩outcome UKB overlap | no UKB-excluded CPD reachable; overlap broken from outcome side in B (Kunkle/Lambert zero-UKB, effect persists) | overlap not the driver; MRlap (full-sumstats) is the next step |
| **4e (D) — REDONE genome-wide** | SlopeHunter index-event correction (lifespan→AD), full sumstats, **2,346 LD-clumped SNPs** (was 10) | b_SH = **−0.944, 95% CI [−1.081, −0.806]** (excl APOE). Corrected at **all three levels**: composite 20-SNP **−0.113 → −0.012 (p=0.73)**; pooled nAChR **−0.129 → +0.005 (p=0.91)**; 15q25 **−0.105 → +0.019 (p=0.70)** | **correction is identified and ABOLISHES the protective effect everywhere** — supersedes the earlier "inconclusive" result |
| **4e (D) coherence** | Does the collider structure quantitatively predict the signal? | predicted induced effect = b_SH × (smk→lifespan): **−0.101** composite / **−0.134** pooled / −0.123 at 15q25 vs **observed −0.113 / −0.129 / −0.105** | the selection pathway accounts for ~**the entire** observed effect at every level (89–117%) |

### Verdict — is the nAChR-locus protective signal separable from selection/collider bias?

> **⚠️ VERDICT REVISED (2026-07-18) after the genome-wide index-event correction.** The earlier
> reading — "collider bias not demonstrated" — rested on Task D being *unidentified* (n=10
> SNPs). Redone properly with full summary statistics and 2,346 LD-clumped SNPs, the correction
> is **tightly identified and it abolishes the protective effect**. The weight of evidence has
> shifted: the signal now looks **largely attributable to survival/competing-risk collider
> bias**. The narrative below is retained for the audit trail, with the revision marked.

**Current reading: the protective nAChR→AD signal is largely or wholly explained by
survival/competing-risk collider bias.** Four independent strands now align — both collider
arms are present (X→S, Y→S), the outcome carries a survival signature (E), the identified
index-event correction removes the entire effect (D), and the collider structure
*quantitatively predicts* the observed magnitude (−0.134 predicted vs −0.129 observed). The
one apparently-contrary result (Task B, persistence in clinical-only AD) tested a **different
channel** — proxy ascertainment — and does not speak to survival selection, which operates in
clinical cohorts too (cases must survive to diagnosis age).

Reading the tests by what they can and cannot establish:

- **Proxy-ascertainment collider — EXCLUDED (Task B, the cleanest test).** If proxy-majority
  Bellenguez were manufacturing the protection, de-selecting the outcome should attenuate it
  toward null. It does the opposite: clinical-only, zero-UKB **Kunkle −0.136 and Lambert
  −0.123 are as protective or more** than proxy-majority Bellenguez −0.105 (difference p=0.74;
  pooled p=0.85). The point estimates move the wrong way for a proxy artefact.
- **Both collider arms are now confirmed present — the structural precondition is fully met,
  but it is *necessary, not sufficient* (Tasks A + Y→S arm + E).** X→S: the 15q25 instruments
  associate with lung cancer (p=10⁻⁵¹) / lifespan / COPD (though 15q25 *is* the lung-cancer
  locus, so this is near-tautological). Y→S: AD liability **shortens lifespan**, and critically
  the **non-APOE** component — the part relevant to the receptor loci, since APOE is absent from
  the nAChR instruments — is itself significant (b=+0.020, p=0.003), albeit modest. Task E adds
  that the outcome shows a survival signature (lung cancer & CAD spuriously protective against
  AD). So a survival/competing-risk collider is a **fully-formed structural possibility**, not
  just a hand-wave — yet "both arms present" establishes only that the collider *can* operate,
  not that it *does*: the induced-bias magnitude depends on joint selection intensity, which
  Tasks B and D probe.
- **Index-event correction — NOW THE STRONGEST EVIDENCE, and it points to a collider
  (Task D, redone genome-wide).** With full summary statistics (Bellenguez AD × Pilling
  lifespan, merged on rsID, allele-aligned, plink2-clumped against 1000G EUR) the fitting set
  goes from **10 → 2,346 independent SNPs** and the slope becomes tightly identified:
  **b_SH = −0.944, 95% CI [−1.081, −0.806]**. Applying it collapses the effect at **every
  level of the decomposition** — composite 20-SNP **−0.113 → −0.012 (p=0.73)**, pooled nAChR
  **−0.129 → +0.005 (p=0.91)**, 15q25 **−0.105 → +0.019 (p=0.70)**. This matters: it is not
  only the receptor *localisation* that dissolves but the original composite signal the whole
  project set out to decompose. Robust
  across variants: distance-pruned replicate (1,178 SNPs) gives −0.876; excluding APOE
  *strengthens* the slope (−0.863 → −0.944); excluding 15q25 changes nothing (−0.9436 →
  −0.9444). **My earlier "APOE/15q25 contamination" diagnosis is retracted** — the n=10 failure
  was purely a power problem, not contamination.
- **Quantitative coherence — the collider predicts the effect size.** Parameterised entirely
  independently (X→S arm from Task A × the SlopeHunter slope), the *predicted* selection-induced
  effect is **−0.134** (pooled) and −0.123 (15q25), against **observed −0.129** and −0.105. The
  survival pathway does not merely attenuate the signal; it accounts for essentially all of it.

The pre-registered "selection excluded" standard (A-equivalence **AND** B-persistence) is
**not met** — A fails on equivalence — and Task D now supplies positive evidence in the other
direction. Task B remains valid but narrower than it first appeared: it excludes the **proxy**
channel only.

### The one caveat that keeps this short of proof

The index-event correction removes the component of the AD effect that is **collinear with the
mortality axis**. It cannot, by construction, distinguish:

1. a **collider-induced** association (what we want to remove), from
2. a **genuine** smoking→AD effect that happens to run along the same axis.

Because smoking strongly affects mortality, *any* true smoking→AD effect would be substantially
collinear with lifespan and would be attenuated by this correction. SlopeHunter's identifying
assumption is that the "incidence-only" SNP cluster reveals a *selection* slope; if lifespan and
AD share genuine biology (they do — vascular, inflammatory, APOE), part of b_SH is biological
and the correction over-removes. So this is strong evidence, not proof.

**Bottom line:** the nAChR-locus protective signal is **most plausibly a survival/competing-risk
collider artefact**. An identified index-event correction eliminates it (−0.129 → +0.005), both
collider arms are present, the outcome carries a survival signature, and the collider structure
quantitatively reproduces the observed effect size. Proxy-ascertainment bias is *not* the
mechanism (Task B) — survival selection is. This does not reach proof, because the correction
cannot separate genuine same-axis pleiotropy from collider-induced association. The decisive
remaining test is unchanged and is now higher priority: an **incident / younger-onset AD**
outcome that does not condition on survival to diagnosis age. **Until that is run, the protective
smoking/nAChR→AD effect should not be treated as a causal, druggable signal.**

## Recommended next steps

> **Priority reordered 2026-07-18.** The brain-QTL work was motivated by localising a
> protective effect that the genome-wide index-event correction now shows is most plausibly
> selection-generated. It is therefore **gated**, not cancelled (APPROACH.md §14).

1. **Incident / younger-onset AD outcome — THE GATE.** The one test that breaks the
   survival-to-diagnosis collider, because such cohorts do not condition on living long enough
   to be ascertained. Not reachable in OpenGWAS; needs individual-level UKBB or a published
   EOAD / incident-dementia GWAS. **Two branches:** effect stays null ⇒ the decomposition
   question closes and this becomes a methodological/negative result; effect survives ⇒ the
   localisation question is live again and steps 2–4 resume with a far stronger rationale.
2. ⛔ **GATED on (1) — Brain QTL / pQTL for the receptor panel.** MetaBrain, GTEx v8 brain,
   ROSMAP eQTL, UKB-PPP / deCODE pQTL in `04_direct_effect.qmd`; CHRNA4/CHRNB2 (α4β2) and
   CHRNA7 (α7) the priority subunits. Design preserved; **do not start before (1)** — it would
   most likely buy an expensive null, and a hit would be uninterpretable (indistinguishable
   from noise fitted to an artefact).
3. ⛔ **GATED on (1)** — CHRNA7 with CHRFAM7A sensitivity (§9); per-locus direction-of-effect
   to convert a surviving target into an agonist-vs-antagonist call (APPROACH.md §7).
4. **Optional corroboration, not gated:** UKB-excluded CPD exposure + MRlap; full CAUSE with
   genome-wide sumstats (now cached under `data/cache/sumstats/`).
5. ~~Re-evaluate the proxy-AD (Wightman) comparison on a harmonised effect scale.~~ **Done in
   L4c** (and the "Wightman" ID was mislabelled). The decisive remaining selection-bias probe is
   an **incident / younger-onset AD** outcome — not reachable in OpenGWAS — to break the
   survival-to-diagnosis collider that L4b/L4c cannot.

## Honesty flags (carried from APPROACH.md §12)

- The novelty here is the **decomposition**, not a generic cis scan — and the decomposition
  is currently incomplete for lack of brain QTL.
- eQTL ≠ channel function; a real functional effect can be invisible to expression
  instruments.
- MVMR with shared cis instruments is prone to conditional weak-instrument bias; conditional
  F could not be estimated for the single-SNP CHRNB2 instrument and is reported as NA rather
  than over-interpreted.
- **The prior L4 selection check tested only CHRNB2, the ≈null contributor — not the
  effect-driving 15q25 locus.** Corrected in L4b: the driving instruments are *not*
  selection-clean (15q25 → lung cancer p=10⁻⁵¹), so the earlier "survival axis amputated"
  claim held only for the one clean locus and did not generalise.
- **Two config accession IDs were mislabelled** (`ieu-a-1001` = Okbay EA not lifespan;
  `ieu-b-5067` = Woolf AD 954 cases not Wightman). The earlier SUMMARY "proxy (Wightman)"
  ascertainment comparison therefore used the wrong dataset — superseded by the verified
  ascertainment gradient in L4c.
- **The collider verdict changed twice — the full trail is recorded deliberately.** (i) The
  original n=10 index-event run said "effect abolished"; (ii) we retracted that, diagnosing
  APOE/circularity contamination; (iii) the genome-wide redo (2,346 SNPs) shows the correction
  is **identified** and *does* abolish the effect, and that the contamination diagnosis was
  **wrong** — excluding APOE strengthens the slope, excluding 15q25 changes nothing. So (i) was
  right for the wrong reason (accurate point estimate, no precision), and (ii) was a
  misdiagnosis of a pure power problem. Current verdict: the signal is **most plausibly a
  survival-collider artefact**, short of proof only because the correction cannot separate
  genuine same-axis pleiotropy from collider-induced association.
- **Task B is valid but narrower than it first read.** Persistence in clinical-only AD excludes
  the **proxy-ascertainment** channel; it does *not* address survival-to-diagnosis selection,
  which operates in Kunkle/Lambert too. It should not be cited as evidence the effect is
  collider-free.
- **Outcome-orientation was verified, not assumed (Y→S arm).** The Pilling parental-longevity
  outcome is a Martingale-residual scale where **higher = shorter life**, so a *positive* MR
  slope means life-shortening. An auto-generated label initially read this backwards
  ("no shortening"); it was corrected after confirming the orientation against known anchors in
  the same dataset (ApoE4 rs429358 β=+0.057, CHRNA5 rs16969968 β=+0.025 — both life-shortening
  alleles are positive). All signs on the lifespan outcome (including Task A) follow this
  convention.
- **Task D is now identified and is the strongest result in the battery** (genome-wide redo:
  full Bellenguez × Pilling sumstats, 148,162 merged SNPs → 2,346 LD-clumped; b_SH = −0.944,
  95% CI [−1.081, −0.806]; corrected pooled −0.129 → +0.005). It **supersedes** the earlier
  n=10 "unidentified/inconclusive" result. The residual limitation is conceptual, not
  statistical: the correction cannot distinguish a collider-induced association from a genuine
  effect collinear with the mortality axis.
- **The remaining uncertainty is now concentrated in one place:** whether b_SH is a pure
  *selection* slope or partly *biological* (lifespan and AD share vascular/inflammatory/APOE
  biology). If partly biological, the correction over-removes and the true effect is somewhere
  between +0.005 and −0.129. An incident/younger-onset AD outcome resolves this; nothing in the
  current data does.
- **Tasks C and D are data-limited by the OpenGWAS API:** a UKB-excluded CPD GWAS, full
  genome-wide summary statistics for MRlap, and an incident/younger-onset AD outcome were not
  reachable and are reported as NA gaps rather than substituted or fabricated.
- The GWAX (Schwartzentruber family-history) outcome in L4c is on a different scale and is
  rescaled ×2 (approximate); the decisive L4c contrast (Bellenguez vs Kunkle/Lambert) is
  scale-clean and does not depend on that rescale.
