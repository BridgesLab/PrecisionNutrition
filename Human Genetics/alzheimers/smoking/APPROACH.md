# APPROACH.md — Decomposing the Smoking→Alzheimer's MR Signal into Druggable nAChR Targets

> Drop this in the project root for Claude Code. It is the analysis contract: scope,
> conventions, pipeline, decision rules, and output spec. Implement each layer as a
> separate Quarto `.qmd`. Prefer frequentist estimators by default; Bayesian variants
> are listed where they add value.

---

## 1. Objective

A protective Mendelian randomization (MR) effect of cigarette smoking on Alzheimer's
disease (AD) is observed. The goal is **not** to confirm that smoking is protective —
it is to **decompose that signal** and isolate the portion that is (a) mediated by
nicotinic acetylcholine receptor (nAChR) *function* and therefore (b) **druggable**, then
nominate candidate therapeutics with a coherent direction of effect.

## 2. Core thesis (do not lose this framing)

The smoking instrument rides two distinct causal axes:

- **Behavioral / combustion axis** — genotype → smoking behavior → nicotine + combustion
  toxins → vascular / inflammatory / mortality effects → AD. This axis carries the
  horizontal pleiotropy **and** the survival/selection collider. **Not druggable** (a
  drug does not make a person smoke).
- **Pharmacodynamic receptor axis** — genotype → nAChR function in brain → direct CNS
  effect on AD, independent of smoking behavior. **This is the only druggable axis**;
  it is what an agonist / PAM recapitulates.

**The decomposition target is the receptor-axis direct effect.** Isolating it also
removes the survival-bias collider, because that collider lives entirely on the
behavioral axis. The druggable estimate and the bias-free estimate are the same quantity.

## 3. Conventions

- **Language/stack:** R + tidyverse/dplyr; analyses authored as Quarto `.qmd` with
  `{r}` chunks. Keep code portable and reusable; no absolute paths (use `here::here()`).
- **Estimation default:** frequentist. Provide Bayesian alternatives where flagged.
- **Reproducibility:** set a seed; pin package versions in `renv.lock`; cache GWAS/QTL
  pulls under `data/cache/`.
- **OpenGWAS auth:** `ieugwasr` uses a JWT (`ieugwasr::get_opengwas_jwt()` /
  `OPENGWAS_JWT` env var). Fail loudly if absent.
- **Significance/thresholds:**
  - Univariable instrument strength: F > 10.
  - MVMR conditional instrument strength: Sanderson–Windmeijer conditional F > 10.
  - Colocalization: `PP.H4 ≥ 0.8` (report H3 vs H4).
  - Steiger directionality: TRUE and p < 0.05.
  - Multiplicity: Bonferroni across the gene panel (report nominal + adjusted).
- **Outputs:** every layer writes a tidy results table to `results/` and a decisions
  line to `results/decisions_log.md`. Plots to `figures/`.

## 4. Key packages

```r
# MR core
library(TwoSampleMR)   # harmonise, IVW/Egger/median, Steiger
library(ieugwasr)      # OpenGWAS pulls, ld_clump
library(MVMR)          # multivariable MR, conditional F (Sanderson–Windmeijer)
library(mrclust)       # MR-Clust mechanism discovery
library(cause)         # Bayesian correlated-pleiotropy check
library(coloc)         # coloc.abf + coloc.susie (runsusie)
library(susieR)        # fine-mapping for SuSiE-coloc
library(tidyverse); library(here)
```

## 5. Data inputs (fill in IDs in `config.yml`)

| Role | Source | Notes |
|---|---|---|
| AD outcome (primary) | Bellenguez 2022 (direct-case) | largest case-ascertained |
| AD outcome (replication) | Wightman 2021; FinnGen | proxy/GWAX — different selection structure (use to probe survival bias) |
| Smoking exposure | GSCAN / Liu 2019 (CigPerDay, SmkInit) | behavioral phenotype |
| Receptor cis-instruments | MetaBrain, ROSMAP, GTEx v8 brain (eQTL); eQTLGen (blood) | brain-tissue preferred |
| Receptor cis-pQTL | UKB-PPP, deCODE | better specificity than eQTL where available |
| Mortality / lifespan | parental lifespan (Timmers), cause-specific (lung cancer, CAD, COPD) | for orthogonality check |
| Pleiotropy covariates | BMI, LDL/total cholesterol, educational attainment | MVMR adjustment set |

**Gene panel (mechanism bins):**

```r
nachr_genes <- c("CHRNA5","CHRNA3","CHRNB4",   # 15q25 cluster
                 "CHRNA4","CHRNB2",            # high-affinity α4β2
                 "CHRNA7",                     # α7 (handle CHRFAM7A — see §9)
                 "CHRNA6","CHRNB3")            # 8p11
mechanism_bins <- list(
  nAChR_pharmacodynamic = nachr_genes,
  nicotine_metabolism   = c("CYP2A6","CYP2B6"),
  reward_behavioral     = c("DRD2","DBH"),
  other_pleiotropic     = NULL   # everything else
)
```

---

## 6. Pipeline

### Layer 0 — Establish the composite signal
Reproduce smoking→AD MR with the full instrument set (IVW + Egger + weighted median +
MR-PRESSO). Record sign, magnitude, heterogeneity. This is the baseline to be decomposed,
not the result.

### Layer 1 — Discovery: is protection mechanistically clustered?
Run **MR-Clust** on the smoking→AD instruments. A real receptor effect should concentrate
in a cluster enriched for nAChR loci, not smear across all SNPs.

```r
cl <- mrclust::mr_clust_em(
  theta    = h$beta.outcome / h$beta.exposure,
  theta_se = abs(h$se.outcome / h$beta.exposure),
  bx = h$beta.exposure, by = h$beta.outcome,
  bxse = h$se.exposure, byse = h$se.outcome
)
# Annotate SNPs → mechanism bin (nearest gene), join cluster assignment,
# then hypergeometric test: is the protective cluster enriched for nAChR loci?
```
**Decision:** if a protective cluster is enriched for `nAChR_pharmacodynamic` → proceed.
If protection sits in metabolism/behavioral/other bins → the signal is dose/combustion,
not a receptor drug; document and stop or reframe.

### Layer 2 — Mechanism-stratified MR
Partition the instrument by bin; run MR per bin. The nAChR-only estimate is the candidate
druggable signal; other bins are non-druggable comparators.

```r
strat <- h_binned |>
  group_by(mechanism) |>
  group_modify(~ mr(.x, method_list = c("mr_ivw","mr_egger_regression",
                                        "mr_weighted_median")) |>
                 select(method, b, se, pval)) |>
  ungroup()
```
**Decision:** effect localizing to `nAChR_pharmacodynamic` supports a receptor mechanism.

### Layer 3 — Formal direct effect (the druggable estimate)

> ⛔ **GATED as of 2026-07-18 — do not start the brain-QTL/pQTL expansion yet.** Layer 4e's
> genome-wide index-event correction removes the protective effect at **every** level
> (composite −0.113 → −0.012; pooled nAChR −0.129 → +0.005; 15q25 −0.105 → +0.019), and the
> selection structure quantitatively reproduces 89–117% of the observed signal. Layer 3 asks
> *"which subunit carries the protective effect?"* — a question that presupposes there is an
> effect to attribute. Running it now would most likely buy an expensive null, and a *hit*
> would be uninterpretable (indistinguishable from noise fitted to an artefact).
>
> **Gate:** run the **incident / younger-onset AD** test first (§14). If the effect survives
> there, Layer 3 becomes live again with a real target and a stronger rationale than it ever
> had. If it stays null, the decomposition question closes and the brain QTL is never spent.
> The design below is preserved unchanged for that first branch.

For each flagged gene, **switch the exposure from the behavioral phenotype to the gene's
cis-eQTL/pQTL**, then partition via MVMR / network MR with smoking heaviness as mediator.
The **direct** (behavior-independent) path is what a drug reproduces.

```r
# Exposures: X1 = receptor cis-expression/protein, X2 = smoking heaviness
# Instruments: receptor cis-SNPs (+ genome-wide smoking SNPs at non-receptor loci)
mv  <- MVMR::format_mvmr(BXGs = cbind(bx_receptor, bx_smk),
                         BYG  = by_AD,
                         seBXGs = cbind(sx_receptor, sx_smk),
                         seBYG  = sy_AD, RSID = rsid)
sf  <- MVMR::strength_mvmr(mv, gen_cov = 0)   # conditional F — the real bottleneck
est <- MVMR::ivw_mvmr(mv)                      # direct effect of receptor | smoking
```
Then **confirm at each locus** with cis-MR (Wald/IVW on cis-SNPs) **+ colocalization**:

```r
co <- coloc::coloc.abf(dataset1 = eqtl_region, dataset2 = ad_region)
# If LD is complex, use coloc.susie(runsusie(...)) to allow multiple causal signals
```
**Decision (target nomination):** keep a gene only if direct-effect MVMR is non-null in
the protective direction **and** `PP.H4 ≥ 0.8`. Conditional F < 10 ⇒ flag weak-instrument
bias, do not over-interpret.

---

## 7. Direction-of-effect rule (the trap — apply per locus)

Sign relationships are locus-specific and counterintuitive. Example: `rs16969968`
(CHRNA5) is a **reduced-function** variant that **increases** cigarettes/day, so "more
smoking" can map to *lower* receptor function. **Never map the behavioral sign to a drug
direction.** For each locus, trace:

1. variant allele → effect on receptor expression/function (sign A, from cis-eQTL/pQTL);
2. same allele → effect on AD (sign B).

Then: higher receptor function → lower AD ⇒ **agonist / PAM** indicated; higher function
→ higher AD ⇒ **antagonist** indicated. Record sign A, sign B, and the implied modality
per gene in the targets table.

## 8. Robustness layer

- **Survival/selection (collider):** look up each receptor instrument against lifespan +
  cause-specific mortality (lung cancer, CAD, COPD). Near-null ⇒ behavioral axis
  successfully amputated. Also compare direct-case (Bellenguez) vs proxy (Wightman)
  outcomes; divergence by ascertainment age is a selection signature.
- **Correlated pleiotropy:** `cause::cause()` on the behavioral instruments.
- **Measured pleiotropy:** MVMR conditioning on BMI / lipids / educational attainment.
- **Directionality:** Steiger filtering at every cis-locus.
- **MVMR Q / heterogeneity:** `MVMR::pleiotropy_mvmr`.

## 9. CHRNA7 / CHRFAM7A special handling

`CHRFAM7A` is a human-specific partial duplication of `CHRNA7` acting as a dominant-
negative regulator; it is a CNV adjacent to the locus and contaminates both eQTL signal
and biology. For any CHRNA7 result: (i) check whether instruments tag the duplication;
(ii) run coloc with and without the duplicated region; (iii) report α7 conclusions as
sensitivity-qualified.

## 10. Target → drug mapping (post-nomination)

| Gene | Modality (if protective = ↑function) | Example agents | Human context |
|---|---|---|---|
| CHRNA7 | α7 agonist / PAM | encenicline, GTS-21, AVL-3288 (PAM) | multiple AD trials failed/discontinued — MR adjudicates target vs molecule |
| CHRNA4 / CHRNB2 | α4β2 agonist / partial agonist | nicotine (transdermal), varenicline | MCI cognition signal; 2-yr MIND trial completed, results pending |
| CHRNA5 | accessory subunit | (no selective ligand) | strongest instrument, weak druggability |
| CHRNA3 / CHRNB4 | α3β4 | varenicline, cytisine (non-selective) | peripheral; little AD data |

Map only genes that survive Layer 3 + §7 direction logic. Flag repurposing-ready agents
(varenicline, nicotine, galantamine) separately from discontinued chemical matter.

## 11. Deliverables Claude Code should produce

1. `01_baseline.qmd` … `05_robustness.qmd` (one per layer) + rendered HTML.
2. `results/targets.tsv` — gene, direct-effect estimate (CI), coloc PP.H4, sign A, sign B,
   implied modality, conditional F, mortality-orthogonality flag, CHRFAM7A flag.
3. Forest plot of mechanism-stratified + per-gene direct effects.
4. `results/decisions_log.md` — one line per decision rule applied, with the value met.

## 12. Honesty flags (state these in the writeup)

- Agnostic druggable-genome MR scans likely already tested naive cis-expression→AD for
  these genes and did **not** prioritize them; the novelty is the **decomposition**, not
  a generic cis scan. Position against those scans explicitly.
- eQTL instruments capture *expression*, not *channel function*; a real functional effect
  can be invisible to an eQTL. Note this limitation.
- MVMR with shared cis instruments is prone to conditional weak-instrument bias — report
  conditional F honestly rather than over-claiming a clean direct effect.

---

## 13. Selection / collider-bias work — pre-registration (added 2026-06-16)

Added to address the open inferential threat that the protective nAChR→AD signal is a
survival/selection collider artefact, not a causal effect. Implemented as Layers 4b–4f
(`04b`–`04f` qmds). The pre-registration below is fixed **before** any Task-A/B test is run.

### 13.1 Scale
All effects on **per-SD cigarettes/day → log-odds outcome**. Selection-axis outcomes
(Task A) are read on their native per-allele scale and the MR slope is interpreted as
log-OR (or SD-units for the quantitative lifespan / FEV1-FVC outcomes) per SD CPD; the
quantitative outcomes are therefore on a *different unit* than the AD log-OR and the
equivalence bound is applied on each outcome's own scale (see 13.3).

### 13.2 Pre-registered equivalence bounds (δ) — TOST
The reference effect to beat is the AD signal itself: pooled nAChR L6 IVW **b = −0.129**
(15q25 b = −0.105). A selection pathway large enough to *manufacture* that AD effect would
have to move the selection-axis outcome by a comparable amount. We therefore pre-register:

- **Primary δ = 0.02** (log-OR / SD-unit per SD CPD) — ~15% of the AD effect magnitude;
  an instrument–selection association below this is too small to plausibly drive the −0.129
  AD estimate.
- **Sensitivity δ = 0.05** — ~40% of the AD effect; a more permissive bound.

A locus×outcome is declared **equivalence-null** only if the **95% CI ⊂ ±δ** (equivalently
TOST p < 0.05 against both bounds), *not* merely p > 0.05 for the point estimate. A wide CI
that merely includes 0 is reported as **"inconclusive"**, never as "clean".

### 13.3 Acceptance rules
- **Task A (selection-clean):** a locus is selection-clean only if equivalence (CI ⊂ ±δ)
  holds for **parental lifespan AND lung cancer AND COPD/FEV1-FVC** simultaneously. 15q25 is
  flagged explicitly as equivalence-null vs real-association vs inconclusive.
- **Task B (ascertainment swap):** substantial attenuation toward null in clinical-only
  outcomes (Kunkle/Lambert) vs the proxy-majority Bellenguez ⇒ supports collider/selection
  bias. Persistence at similar magnitude with comparable CIs ⇒ supports a genuine effect.
  Formal test: z on the effect difference; meta-regression of effect on proxy fraction
  across the ≥3-outcome gradient (proxy_frac ≈ 0 Kunkle/Lambert, 0.54 Bellenguez, 1.0
  Schwartzentruber-GWAX).
- The standard for "selection bias excluded" is **A (equivalence) AND B (clinical-only
  persistence) together** — neither a single instrument-on-lifespan null nor B alone suffices.

### 13.4 Verified accessions (composition-audited 2026-06-16; see
`results/dataset_composition_verified.csv`)
Two original config IDs were found MISLABELLED and are **not** used as labelled:
`ieu-a-1001` ("lifespan") is actually Okbay EA; `ieu-b-5067` ("Wightman proxy") is actually
Woolf 2022 AD with 954 cases. Verified replacements live in `config.yml` under
`selection_axis:` and `ad_ascertainment:`. Unreachable targets (all-cause-mortality binary,
UKB-participation GWAS, incident/younger-onset AD) are reported as **NA gaps**, not
substituted.

---

## 14. Outcome and revised priority (2026-07-18)

### 14.1 What the selection work concluded
The index-event correction was redone with **full genome-wide summary statistics** (Bellenguez
GCST90027158 × Pilling GCST006697, merged on rsID, allele-aligned, `plink2 --clump` r²<0.01/1Mb
against 1000G EUR), taking the SlopeHunter fitting set from **10 → 2,346 independent SNPs**. The
slope is tightly identified (**b_SH = −0.944, 95% CI [−1.081, −0.806]**) and the correction
**abolishes the protective effect at every level of the decomposition**:

| Set | uncorrected | corrected | p | % of signal explained by selection |
|---|---|---|---|---|
| canonical 20-SNP composite | −0.113 | −0.012 | 0.73 | 89% |
| pooled nAChR | −0.129 | +0.005 | 0.91 | 104% |
| 15q25 CHRNA5/A3/B4 | −0.105 | +0.019 | 0.70 | 117% |

The collider structure, parameterised independently as `b_SH × (smk→lifespan)`, reproduces the
observed magnitude in each case. Both collider arms are confirmed present (X→S, Y→S) and the
outcome carries a survival signature. **Proxy ascertainment is not the mechanism** (the effect
persists in clinical-only AD, Layer 4c) — **survival/competing-risk selection is.**

### 14.2 The one caveat short of proof
The correction removes whatever is **collinear with the mortality axis** and cannot separate a
collider-induced association from a *genuine* smoking→AD effect running along that same axis.
Since smoking strongly affects mortality, a true effect would also be attenuated. If b_SH is
partly biological (lifespan and AD share vascular/inflammatory/APOE biology), the truth lies
between the corrected and uncorrected values.

### 14.3 Revised priority — the gate
1. **Incident / younger-onset AD outcome** (does not condition on survival to diagnosis age).
   This is now the single decisive test and the gate on everything downstream. Not reachable in
   OpenGWAS; needs individual-level UKBB or a published EOAD/incident-dementia GWAS.
2. **Gated on (1):** Layer 3 brain eQTL/pQTL + colocalization (§6) — preserved but not started.
3. Optional corroboration: UKB-excluded CPD exposure + MRlap; APOE-excluded SlopeHunter
   variants (already run as sensitivity — conclusion unchanged).

**Until (1) is run, the protective smoking/nAChR→AD effect should not be treated as a causal,
druggable signal.**
