# Nicotine Drug-Target MR — Decomposing the Smoking → Alzheimer's Signal

This project implements the analysis contract in [`APPROACH.md`](APPROACH.md): take the
observed **protective** Mendelian-randomization (MR) effect of cigarette smoking on
Alzheimer's disease (AD) and try to **decompose** it into a *druggable* nicotinic
acetylcholine receptor (nAChR) pharmacodynamic axis vs. a *non-druggable*
behavioral/combustion axis. The deliverable is a set of layered Quarto analyses, tidy
results tables, and a target-nomination table.

See [`SUMMARY.md`](SUMMARY.md) for the results and interpretation.

## ⚠️ Headline finding (2026-07-18) — read before using any downstream result

The decomposition **worked**, but the signal it decomposed does not survive selection-bias
correction. Layers **4b–4e** (added after the original 01–06 pipeline) test the
survival/selection collider directly:

- Both collider arms are present — the effect-driving 15q25/nAChR instruments engage the
  mortality axis (lung cancer p=10⁻⁵¹), **and** AD liability shortens lifespan (incl. a
  significant non-APOE component).
- A **genome-wide index-event correction** (SlopeHunter, 2,346 LD-clumped SNPs;
  b_SH = −0.944, 95% CI [−1.081, −0.806]) **abolishes the protective effect at every level**:
  composite 20-SNP −0.113 → −0.012; pooled nAChR −0.129 → +0.005; 15q25 −0.105 → +0.019.
- The collider structure *quantitatively* reproduces **89–117%** of the observed signal.

**Consequence:** the protective smoking/nAChR→AD effect should **not** be treated as a causal,
druggable signal, and the Layer 3 brain-QTL expansion is **gated** (APPROACH.md §14) pending an
**incident / younger-onset AD** replication. Proxy ascertainment is *not* the mechanism
(Layer 4c) — survival selection is.

## How to run

Requirements:

- R ≥ 4.4 with: `TwoSampleMR`, `ieugwasr`, `MVMR`, `MendelianRandomization`, `mrclust`,
  `coloc`, `susieR`, `cause`, `MRPRESSO`, `mr.raps`, `tidyverse`, `here`, `yaml`,
  `GenomicRanges`, `TxDb.Hsapiens.UCSC.hg19.knownGene`, `org.Hs.eg.db`, `ggrepel`.
  **Layers 4b–4e additionally need:** `data.table`, `SlopeHunter`, `MRlap`, `TOSTER`
  (`SlopeHunter`/`MRlap` are GitHub-only: `remotes::install_github(c("Osmahmoud/SlopeHunter","n-mounier/MRlap"))`).
- A **valid OpenGWAS JWT** in `~/.Renviron` as `OPENGWAS_JWT='...'` (tokens expire — renew
  at <https://api.opengwas.io> if you get HTTP 401).
- **`plink2`** on `PATH` — required for the genome-wide LD clumping in Layer 4e.
- Quarto (e.g. the binary bundled with RStudio at
  `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`).

Render a layer:

```sh
quarto render 01_baseline.qmd --to html
```

Layers 01 → 06 are designed to run **in order**; each reads the previous layer's outputs from
`results/`. Layers **4b–4e** are the selection/collider battery and read the L6 cis-instrument
caches plus (for 4e) the large external inputs below.

### External data (Layer 4e) — provenance and how to re-fetch

Layer 4e needs **full genome-wide summary statistics**, which the OpenGWAS API cannot serve (it
returns clumped tophits and regional queries only — that limitation is exactly why the first
attempt at Task D was underpowered at n=10). One script fetches everything (~3.5 GB, idempotent):

```sh
bash scripts/fetch_sumstats.sh
```

| File | Source | Role |
|---|---|---|
| `data/cache/sumstats/AD_bellenguez_harmonised.tsv.gz` | GWAS Catalog **GCST90027158** (Bellenguez 2022, harmonised, GRCh38) | prognosis trait (AD) |
| `data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz` | GWAS Catalog **GCST006697** (Pilling 2017 parental longevity, harmonised) | incidence / selection axis |
| `data/cache/EUR.{bed,bim,fam}` | 1000 Genomes EUR panel, MRC-IEU mirror (`1kg.v3.tgz`) | LD reference for `plink2 --clump` |

> **Sign convention — important.** The Pilling outcome is on a **Martingale-residual** scale
> where **higher = SHORTER lifespan**, so a *positive* MR slope means life-shortening. Verified
> against known anchors in that dataset: ApoE4 `rs429358` β=+0.057 and CHRNA5 `rs16969968`
> β=+0.025 (both life-shortening alleles are positive). All lifespan-outcome signs in Layers
> 4b/4e follow this convention.

The script prints the three follow-on commands (merge/align → `plink2 --clump` → SlopeHunter fit).
Once `results/indexevent_*.csv` exist, `04e_overlap_indexevent.qmd` renders **offline** — its
live chunks read those CSVs and need no network.

## File inventory

### Configuration & helpers
| File | Purpose |
|---|---|
| [`config.yml`](config.yml) | OpenGWAS dataset IDs, sample sizes, thresholds (F>10, PP.H4≥0.8, Steiger p<0.05), and the nAChR + comparator gene panel / mechanism bins. Edit IDs here. |
| [`R/helpers.R`](R/helpers.R) | Shared functions: nearest-gene annotation (hg19); **locus-based mechanism binning** (`assign_mechanism_bin` / `panel_gene_coords` — assigns SNPs by panel-gene cis-window membership with a nearest-gene fallback, recording `bin_gene` + `bin_method`); instrument-strength (R²/F); LD-aware cis-MR helpers; `og_retry` OpenGWAS backoff; the `decisions_log` writer. |
| [`APPROACH.md`](APPROACH.md) | The original analysis contract (scope, pipeline, decision rules). Source of truth for the design. |

### Analysis layers (Quarto `.qmd` → rendered `.html`)
| File | Layer | What it does |
|---|---|---|
| [`01_baseline.qmd`](01_baseline.qmd) | **L0 — composite signal** | Reproduces smoking (cigarettes/day, Liu 2019) → AD (Bellenguez 2022) MR: IVW (FE/MRE), RAPS, Egger, weighted median/mode, MR-PRESSO; annotates nearest gene + instrument strength; writes the harmonised+annotated SNP table used by all later layers. |
| [`02_discovery.qmd`](02_discovery.qmd) | **L1 — discovery** | MR-Clust run twice: on the **full** 22-SNP set (confirms CLU/MINDY2 as a distinct risk cluster) and on the **pruned** canonical set (tests for residual substructure); mechanism-bin annotation + hypergeometric nAChR-enrichment test. |
| [`03_stratified.qmd`](03_stratified.qmd) | **L2 — mechanism-stratified MR** | Partitions the **pruned** instrument by **locus-based** mechanism bin (cis-window membership of panel genes, nearest-gene fallback — so 15q25→nAChR, chr19→metabolism) and runs MR within each; forest plot. This is where the protection localises to the nAChR bin (p≈4×10⁻⁴) while "other" goes null. |
| [`04_direct_effect.qmd`](04_direct_effect.qmd) | **L3 — druggable direct effect** | Switches exposure to each receptor gene's **cis-eQTL**; per-gene cis-MR, MVMR (receptor \| smoking) with Sanderson–Windmeijer conditional F, and colocalization (`coloc.abf`) vs AD; builds the target-nomination table with per-locus direction-of-effect → implied drug modality. |
| [`05_robustness.qmd`](05_robustness.qmd) | **L4 — robustness** | Receptor-instrument orthogonality vs lifespan / lung cancer / CAD (survival-collider check); correlated/directional pleiotropy; MVMR conditioning on BMI/LDL/EA; Steiger directionality; direct-case (Bellenguez) vs proxy (Wightman) ascertainment comparison. |
| [`06_cis_drugtarget.qmd`](06_cis_drugtarget.qmd) | **L6 — relaxed cis drug-target MR** | Drug-target-MR convention (HMGCR/PCSK9/NPC1L1 style): relaxes instrument selection within each nAChR cis window (p<1e-5, correlated SNPs r²<0.3) and runs **LD-aware** IVW/Egger (`MendelianRandomization`, variant correlation matrix) for both smoking-cis and cis-eQTL exposures → AD; per-locus + pooled; SNP-gain summary; **15q25 LD-block resolution** (independent-signal count, subunit-attribution limit) + formal CHRNA4-vs-15q25 comparison; strict-vs-relaxed comparison. Correctly captures the 15q25 CHRNA5/A3/B4 cluster the genome-wide set missed. |

### Selection / collider-bias battery (added 2026-06-16 … 2026-07-18)
| File | Layer | What it does |
|---|---|---|
| [`04b_selection_robustness.qmd`](04b_selection_robustness.qmd) | **L4b — the two collider arms** | **Part 1 (X→S):** locus-resolved negative-control battery — the *effect-driving* 15q25/CHRNA4/CHRNB2/CHRNA6B3/pooled instruments vs parental lifespan, lung cancer, COPD(FEV1/FVC), CAD, with **pre-registered TOST equivalence** (δ=0.02 primary, 0.05 sensitivity; APPROACH.md §13). **Part 2 (Y→S):** AD liability (clinical Kunkle instruments) → parental lifespan, with/without APOE. Verdict: both arms present. |
| [`04c_ascertainment_swap.qmd`](04c_ascertainment_swap.qmd) | **L4c — ascertainment swap** | Re-estimates smoking→AD across a **by-proxy gradient** (Kunkle/Lambert proxy=0 → Bellenguez 0.54 → Schwartzentruber GWAX 1.0) on a harmonised per-SD log-OR scale; formal proxy→clinical difference test + meta-regression. Effect persists in clinical-only AD ⇒ **proxy ascertainment excluded**. |
| [`04d_direction_negcontrols.qmd`](04d_direction_negcontrols.qmd) | **L4d — outcome diagnostics** | Direction negative controls through the *same* primary outcome (EA→AD, lung-cancer→AD, CAD→AD). Lung cancer and CAD appear spuriously **protective** ⇒ the outcome carries survival/competing-risk structure. |
| [`04e_overlap_indexevent.qmd`](04e_overlap_indexevent.qmd) | **L4e — overlap + index-event correction** | Task C: exposure∩outcome UKB overlap (no UKB-excluded CPD reachable; overlap broken from the outcome side in L4c). Task D **(PRIMARY, genome-wide)**: SlopeHunter index-event correction on 2,346 LD-clumped SNPs → **abolishes the protective effect at all levels**, plus a quantitative coherence check. The original n=10 run is retained as a non-executing legacy section. |

### Inputs
| Path | Contents |
|---|---|
| `data/Instruments - cig.liu - Sleep.csv` | Pre-clumped cigarettes-per-day instrument set (GSCAN/Liu 2019, `ieu-b-142`), copied from the parent Sleep:AD project. |
| `data/cache/instr_*.rds` | Cached L6 relaxed-cis instrument sets (15q25=31, CHRNA4=4, CHRNB2=2, CHRNA6B3=6, pooled=43) — **fixed** so L4b/L4c test the same instruments that produced the protective estimate. |
| `data/cache/sumstats/`, `data/cache/EUR.*` | Large external inputs for L4e — see "External data" above; fetch with `scripts/fetch_sumstats.sh`. |

### Outputs
| Path | Contents |
|---|---|
| `results/baseline_harmonised_annotated.csv` | Harmonised, Steiger-filtered, gene-annotated, bin-assigned SNP table (L0). |
| `results/baseline_mr_results.csv` | Composite MR estimates (all methods) (L0). |
| `results/baseline_pleiotropy_pruned.csv` | IVW/median with vs without the CLU/MINDY2 pleiotropy outliers (L0b). |
| `results/baseline_harmonised_pruned.csv` | **Canonical** 20-SNP annotated/binned set (CLU/MINDY2 removed); read by L1–L5 (L0b). |
| `data/cig_instruments_pruned.csv` | Pruned smoking instrument file (CLU/MINDY2 removed) used by MVMR/CAUSE/ascertainment (L0b). |
| `results/mrclust_assignments.csv`, `mrclust_enrichment.csv` | MR-Clust on the **pruned** set: cluster membership + nAChR enrichment (L1). |
| `results/mrclust_assignments_full.csv`, `mrclust_enrichment_full.csv` | MR-Clust on the **full** set (shows the CLU/MINDY2 risk cluster) (L1). |
| `results/stratified_mr_results.csv` | MR by mechanism bin (L2). |
| `results/cis_mr_per_gene.csv`, `coloc_per_gene.csv`, `mvmr_direct_effect.csv` | Per-gene cis-MR, colocalization, MVMR direct effects (L3). |
| `results/targets.tsv` | **Target-nomination table**: gene, cis effect (CI), coloc PP.H4, sign A/B, implied modality, conditional F, CHRFAM7A flag, nominated y/n (L3). |
| `results/mortality_orthogonality.csv`, `correlated_pleiotropy.csv`, `mvmr_measured_pleiotropy.csv`, `ascertainment_comparison.csv` | Robustness tables (L4). |
| `results/cis_smoking_ldaware_per_locus.csv`, `cis_smoking_ldaware_pooled.csv` | Relaxed smoking-cis LD-aware MR per nAChR locus + pooled (L6). |
| `results/cis_eqtl_ldaware.csv`, `cis_strict_vs_relaxed.csv`, `cis_snp_gain.csv` | Relaxed cis-eQTL MR, strict-vs-relaxed comparison, per-locus SNP-gain table (L6). |
| `results/decisions_log.md` | One line per decision rule applied, with the value met. Periodically de-duplicated; superseded conclusions kept as ⚠️-marked rows for auditability. |
| `results/dataset_composition_verified.csv` | **Composition audit** of every accession used (N, cases/controls, author, year) — after two config IDs were found mislabelled (L4b). |
| `results/selection_battery.csv` | L4b Task A: each nAChR instrument set × selection-axis outcome, with TOST equivalence verdicts at δ=0.02/0.05. |
| `results/collider_Yarm_AD_to_lifespan.csv` | L4b Part 2: AD liability → parental lifespan, with and without APOE. |
| `results/ascertainment_swap.csv`, `ascertainment_metareg.csv` | L4c: effect across the by-proxy gradient + attenuation/meta-regression tests. |
| `results/direction_negcontrols.csv` | L4d: EA / lung-cancer / CAD → AD selection diagnostics. |
| `results/indexevent_slopehunter_fits.csv`, `indexevent_correction_genomewide.csv` | **L4e genome-wide (primary)**: SlopeHunter fits (full / excl-APOE / excl-APOE+15q25) and corrected estimates for composite, pooled nAChR and 15q25. |
| `results/indexevent_correction.csv` | ⚠️ Legacy n=10 correction — superseded; retained for the audit trail. |
| `figures/` | All plots (`.png` + `.pdf`): baseline forest, MR-Clust scatter, stratified forest, selection-battery forest, ascertainment-swap forest. |

## Open TODO — the one dataset that would settle this

> **Wanted: an incident / younger-onset AD (or incident-dementia) GWAS.**
> This is now **the gate on the whole project** (APPROACH.md §14.3). Every test we can run
> conditions on survival to AD-ascertainment age, so the survival collider cannot be broken from
> the inside — including Layer 4c's clinical-only cohorts (Kunkle/Lambert cases still had to live
> long enough to be diagnosed). An **incident** or **early-onset (EOAD)** outcome does not
> condition on that, so it separates the two explanations cleanly:
>
> - effect **stays null** ⇒ the decomposition question closes; this becomes a methodological /
>   negative result and the brain-QTL work is never spent;
> - effect **survives** ⇒ localisation is live again and Layer 3 resumes with a far stronger
>   rationale than it ever had.
>
> **Status:** not found in OpenGWAS — but the search so far was only a handful of guessed
> accessions, **not** a systematic one. Worth a proper GWAS Catalog full-text search
> (EOAD / age-at-onset / incident dementia) before accepting that this needs individual-level
> UK Biobank. If such a dataset turns up, it slots straight into `04c_ascertainment_swap.qmd`
> as another point on the ascertainment gradient (`config.yml: ad_ascertainment.incident_onset`,
> currently `NA`).

## Notes & honest limitations

- **The protective effect is most plausibly a survival-collider artefact** (see headline). The
  one caveat short of proof: the index-event correction removes whatever is *collinear with the
  mortality axis* and cannot separate collider-induced association from a genuine effect running
  along that same axis. If `b_SH` is partly biological (lifespan and AD share
  vascular/inflammatory/APOE biology) the truth lies between corrected and uncorrected values.
- **The collider verdict changed twice** — n=10 "abolished" → retracted as
  contamination → genome-wide "abolished, identified". The full trail is in SUMMARY.md's
  honesty flags and the ⚠️ rows of `results/decisions_log.md`. The middle step was a
  misdiagnosis of what was purely a power problem.
- **Layer 3 is data-limited *and* gated.** Of the 8 panel genes only **CHRNB2** has a blood
  cis-eQTL in OpenGWAS; the rest are brain-predominant and need **brain QTL** (MetaBrain, GTEx
  v8 brain, ROSMAP) or **pQTL** (UKB-PPP, deCODE). That expansion is designed but **should not
  be started** before the incident-AD test above (APPROACH.md §14) — it would most likely buy an
  expensive null, and a hit would be uninterpretable.
- **Reported gaps, not substituted:** all-cause-mortality (binary), a UKB-participation GWAS,
  a UKB-excluded CPD exposure, and incident/younger-onset AD were all unreachable and are
  returned as `NA` rather than proxied with a different dataset.
- eQTL captures *expression*, not channel *function* — a true functional effect can be
  invisible to an eQTL instrument.
- CHRNA7 results are flagged for the `CHRFAM7A` partial-duplication confound (APPROACH.md §9).
- CAUSE in `05_robustness.qmd` uses an MR-Egger-intercept proxy unless full genome-wide
  (un-clumped) summary statistics are cached under `data/cache/`.
