# Decisions log

(Curated in pipeline order; layers append chronologically when re-run, so this file is
periodically de-duplicated. **Last de-duplicated 2026-07-18.** Superseded conclusions are kept
as single ⚠️-marked rows rather than deleted, so the reasoning trail stays auditable — see
SUMMARY.md "Honesty flags" for the two occasions the collider verdict changed.)

| Layer | Rule | Value | Decision |
|---|---|---|---|
| L0 | Composite IVW sign/magnitude (all 22 SNPs) | b=-0.065, p=0.296 | protective but diluted by pleiotropy outliers |
| L0b | IVW after removing CLU/MINDY2 outliers | b=-0.113, p=0.00034 (n=20) | **adopt pruned 20-SNP set as canonical baseline** for L1–L5 |
| L1 | MR-Clust on full set | CLU/MINDY2 = distinct risk cluster (θ=+1.24) | confirms the two are AD-direct horizontal pleiotropy |
| L1 | MR-Clust on pruned set | 1 protective cluster, no further substructure | no other relevant clusters once outliers removed |
| L2 | Mechanism binning rule | cis-window membership + nearest-gene fallback | biology-based bins: 15q25→nAChR, chr19→metabolism (not literal nearest gene) |
| L2 | nAChR-bin IVW (locus-binned, pruned) | b=-0.149, p=0.0004 (n=4) | **protection localises to nAChR**; "other" bin null (b≈0, p=0.95) |
| L3 | Genes nominated (coloc H4>=0.8 & cis-MR p<0.05) | 0 gene(s) | none — blood-eQTL coverage limits brain receptor genes |
| L4 | ⚠️ SUPERSEDED — receptor-instrument mortality orthogonality | 3/3 near-null, but **CHRNB2 only** (the ≈null contributor) | superseded by L4b: the *effect-driving* 15q25 instruments are NOT selection-clean; the "axis amputated" reading did not generalise |
| L6 | Pooled nAChR relaxed-cis LD-aware IVW | b=-0.129, p=0.00087 (n=43) | significant protective receptor-locus effect (drug-target convention) |
| L6 | SNP gain from relaxed criteria | 15q25 0→31, 8p11 1→6 SNPs | newly instruments 15q25 CHRNA5/A3/B4 + CHRNA6/B3; no off-panel targets (cis-restricted) |
| L6 | 15q25 subunit resolution | 7 independent signals, max\|r\|=0.56 | locus-level effect only — CHRNA5/A3/B4 not separable without subunit brain QTL |
| L6 | CHRNA4 vs 15q25 effect difference | diff=-0.15, z=-0.97, p=0.33 | no evidence CHRNA4 stronger than 15q25 |
| L6 | nAChR-cis-eQTL → AD (CHRNB2, relaxed) | b=-0.011, p=0.79 (n=12) | molecular expression effect null; brain eQTL still required to test function |
| L4b | Composition audit of config IDs | ieu-a-1001=Okbay EA (not lifespan); ieu-b-5067=Woolf AD 954 cases (not Wightman) | both MISLABELLED — replaced with verified accessions |
| L4b | Pre-registered equivalence bound δ | δ=0.02 primary, δ=0.05 sensitivity (log-OR/SD per SD CPD) | locus selection-clean iff 95% CI ⊂ ±δ for lifespan AND lung cancer AND COPD |
| L4b | Verified selection-axis outcomes | parental lifespan ebi-a-GCST006697 (Pilling, UKB); lung cancer ebi-a-GCST004748 (McKay, non-UKB); COPD ebi-a-GCST007431 (Shrine FEV1/FVC); CAD ieu-a-7 (Nikpay) | all-cause-mortality + UKB-participation unreachable → NA gaps |
| L4b | Task A (X→S arm): 15q25 vs selection axis (TOST δ=0.02) | lung cancer b=1.36 p=1e-51; parental lifespan b=0.130 p=1e-21; COPD/FEV1FVC b=-0.084 p=5e-11 — all REAL-assoc, none equiv-null | **15q25 is NOT selection-clean**; effect-driving locus engages the mortality axis hard |
| L4b | Task A: any selection-clean locus? | **0/5 sets** equiv-null on lifespan&lungCa&COPD (15q25, CHRNA4, CHRNB2, CHRNA6B3, pooled — all 0/3; CHRNB2 only "inconclusive") | no nAChR set passes equivalence — **X→S arm present** |
| L4b | Task A Part 2 (Y→S arm): AD liability → parental lifespan | all b=+0.043 p=5.5e-31; **excl-APOE b=+0.020 p=0.0032** (Martingale resid: +ve = SHORTER life, verified via APOE rs429358 β=+0.057 & CHRNA5 rs16969968 β=+0.025) | AD liability **SHORTENS lifespan incl non-APOE** → **Y→S arm present**; with X→S, **BOTH collider arms confirmed** (necessary, not sufficient) |
| L4c | Verified AD ascertainment gradient | Kunkle ieu-b-2 & Lambert ieu-a-297 (proxy=0); Bellenguez GCST90027158 (proxy~0.54); Schwartzentruber GCST90012878 GWAX (proxy~1) | enables proxy→clinical attenuation + meta-regression |
| L4c | Task B: proxy→clinical attenuation | 15q25 clinical -0.130 vs Bellenguez -0.105 (Δ z=0.34, p=0.74); pooled -0.117 vs -0.129 (z=-0.19, p=0.85); CHRNA4 p=0.74; canonical 20-SNP p=0.53 | effect does **NOT** attenuate in clinical-only AD — ⚠️ but this excludes the **proxy-ascertainment channel only**, not survival-to-diagnosis selection |
| L4c | Task B: meta-reg effect~proxy_frac | 15q25 slope=+0.140 p=0.037 (driven by rescaled GWAX point); sign is OPPOSITE to collider prediction | proxy ascertainment is not biasing the estimate toward protection |
| L4d | Task E: outcome selection signature | lung cancer→AD OR=0.93 p=0.0036; CAD→AD OR=0.90 p=1.9e-7 (both spuriously protective); EA→AD +0.106 (risk) | **Bellenguez outcome DOES carry survival/competing-risk selection structure** |
| L4e | Task C: exposure-outcome UKB overlap | CPD ieu-b-142 (GSCAN, incl UKB) × Bellenguez controls (~401k, largely UKB) = substantial overlap; no UKB-excluded CPD reachable in OpenGWAS (ieu-b-142 & ieu-b-25 both incl UKB) | overlap broken from OUTCOME side in Task B (Kunkle/Lambert zero-UKB, effect persists); MRlap needs genome-wide sumstats |
| L4e | ⚠️ SUPERSEDED — Task D on OpenGWAS tophits | b_SH=-0.783 (n=10); sensitivity dropping APOE+15q25 gave -0.758 with 95% CI [-1.79,+0.27] **including 0** | read at the time as "unidentified / inconclusive, contamination-driven" — **both that reading and its contamination diagnosis are superseded by L4e-GW** |
| **L4e-GW** | Task D REDONE genome-wide | full sumstats (Bellenguez GCST90027158 × Pilling GCST006697) merged on rsID, allele-aligned, palindromic dropped → 148,162 SNPs; `plink2 --clump` r²<0.01/1Mb vs 1000G EUR → **2,346 independent SNPs** (was 10) | correction is now **IDENTIFIED** |
| **L4e-GW** | SlopeHunter slope (LD-clumped, excl APOE) | **b_SH=-0.944, SE=0.070, 95% CI [-1.081,-0.806]** (n_fit=2,250); incl APOE -0.863; excl APOE+15q25 -0.944; distance-pruned replicate -0.876 | slope tightly estimated; CI excludes 0 in every variant |
| **L4e-GW** | Earlier "APOE/15q25 contamination" diagnosis | excluding APOE **strengthens** the slope (-0.863→-0.944); excluding 15q25 changes it by <0.001 | **RETRACTED** — the n=10 failure was pure **power**, not contamination |
| **L4e-GW** | Index-event-corrected effect, all levels | canonical 20-SNP **-0.113→-0.012** (p=0.73); pooled nAChR **-0.129→+0.005** (p=0.91); 15q25 **-0.105→+0.019** (p=0.70) | **correction ABOLISHES the protective effect at every level** of the decomposition |
| **L4e-GW** | Canonical 20-SNP X→S arm | smk→lifespan b=+0.107, se=0.0148, p=5.3e-13 (n=20; local Pilling sumstats, allele-aligned) | composite instruments engage the mortality axis like the nAChR sets |
| **L4e-GW** | Quantitative coherence of the collider | predicted induced = b_SH × (smk→lifespan) = -0.101 composite / -0.134 pooled / -0.123 15q25 vs observed -0.113 / -0.129 / -0.105 | selection pathway explains **89–117%** of the observed signal at every level |
| **L4e-GW** | Overall verdict | both collider arms present + outcome survival signature + identified correction abolishing the effect + quantitative coherence | protective smoking/nAChR→AD signal is **most plausibly a survival/competing-risk collider artefact**; Layer 3 brain-QTL work **GATED** on an incident/younger-onset AD test (APPROACH.md §14) |
