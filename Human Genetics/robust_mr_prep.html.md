---
title: "Robust MR — Data Preparation and QC"
author: "Dave Bridges"
date: today
editor: source
format:
  html:
    toc: true
    toc-location: right
    keep-md: true
    embed-resources: true    # self-contained: no _files/ dir needed to view it
    code-fold: true
    code-summary: "Show the code"
    fig-path: "figures-robust-mr/"
theme: journal
execute:
  echo: true
  warning: false
---


::: {.cell}

:::


## Purpose

Prepare genome-wide summary statistics for the two overlap-robust MR methods used
in this project: **CAUSE** (Morrison et al. 2020) and **MR-APSS** (Hu et al. 2022).

This is a distinct arm from the existing analyses. `drug_target_mr_bmd.qmd` uses
**cis** instruments at HMGCR/PCSK9/NPC1L1 and asks a drug-target question;
`mr-downstream-analyses.qmd` uses genome-wide clumped instruments with IVW/Egger.
CAUSE and MR-APSS are **genome-wide polygenic** methods — they cannot be applied
to a handful of cis SNPs, and they do not replace the cis analyses. What they add
is a principled treatment of three problems the conventional estimates cannot
handle:

| Problem | Why it bites here | CAUSE | MR-APSS |
|---|---|---|---|
| **Sample overlap** | The LDL exposure and Morris 2019 eBMD are both UK Biobank — effectively one-sample MR | Estimates $\rho$, the null correlation of test statistics, genome-wide | Estimates a 2×2 $C$ matrix from bivariate LDSC intercepts |
| **Correlated horizontal pleiotropy** | Lipids share a large polygenic background with adiposity and with bone | Explicit shared-factor $U$ with parameters $\eta, q$ | Background model ($\Omega$) absorbs it before causal inference |
| **Weak instrument bias** | Relaxing the threshold to gain power normally imports winner's curse | Uses $p<10^{-3}$ variants, modelled rather than filtered | $p<5\times10^{-5}$ **with** an explicit selection-bias correction |

::: callout-important
## Why MRAID is not in this arm

MRAID (Yuan et al., *Sci Adv* 2022) was considered and deliberately **excluded**.
It models correlated and uncorrelated pleiotropy well, but it is a strict
*two-sample* method: its likelihood assumes the exposure and outcome GWAS come
from independent samples. Applying it to UK Biobank LDL against UK Biobank eBMD
would violate the assumption it is least robust to, and its estimate would be
pulled toward the observational association.

MRAID becomes appropriate once we have a genuinely non-overlapping design —
e.g. **UK Biobank exposure → MGI / BioVU outcome**, which is the same
architecture already used for the calcium outcome in `mr-tc-calcium.qmd`.
Recorded as a planned extension, not a gap.
:::

## Study design — two exposure arms

The overlap correction is only credible if we can show it does something. So both
methods are run twice:

- **Overlapping arm** — GLGC 2021 (Graham, EUR, N ≈ 1.32M) → Morris 2019 eBMD.
  Maximum power; UK Biobank is on both sides.
- **Overlap-free arm** — GLGC 2013 (Willer, N ≈ 173k, entirely pre-UK Biobank;
  already in `raw_data/jointGwasMc_LDL.txt.gz`) → Morris 2019 eBMD.
  Lower power, but no overlap by construction.

If the estimates agree after correction and disagree before it, the correction is
doing real work. MR-APSS gives us a bonus here: $C_{12}$, the off-diagonal of the
$C$ matrix, is a *direct empirical measurement* of sample overlap, so we do not
have to take anyone's word for which cohorts went into which meta-analysis.

::: callout-warning
`ANALYSIS.md` §3 currently states that the GLGC total-cholesterol exposure
(`ebi-a-GCST90025953`) is "non-overlapping with UKB". **This needs checking** —
UK Biobank contributes to the Graham 2021 European meta-analysis. Rather than
resolving it from the literature, this pipeline measures it: see the $C_{12}$
estimate in `robust_mr_apss.qmd`.
:::

## Dependencies


::: {.cell}

```{.r .cell-code}
deps <- check_robust_mr_deps(install = FALSE)
deps |> kable(caption = "Package availability for the robust MR arm")
```

::: {.cell-output-display}


Table: Package availability for the robust MR arm

|package     |source               |installed |
|:-----------|:--------------------|:---------|
|tidyverse   |CRAN                 |TRUE      |
|data.table  |CRAN                 |TRUE      |
|yaml        |CRAN                 |TRUE      |
|here        |CRAN                 |TRUE      |
|R.utils     |CRAN                 |TRUE      |
|ieugwasr    |CRAN                 |TRUE      |
|TwoSampleMR |CRAN                 |TRUE      |
|knitr       |CRAN                 |TRUE      |
|kableExtra  |CRAN                 |TRUE      |
|cause       |jean997/cause        |TRUE      |
|MRAPSS      |YangLabHKUST/MR-APSS |TRUE      |
|mixsqp      |stephenslab/mixsqp   |TRUE      |
|ashr        |stephens999/ashr     |TRUE      |


:::
:::


Install the missing ones with:

```r
check_robust_mr_deps(install = TRUE)
# or, explicitly:
# remotes::install_github(c("jean997/cause", "YangLabHKUST/MR-APSS"))
```

## Step 1 — inspect the raw headers before trusting the config

Column names in the GLGC and GEFOS releases have changed between versions. A
mis-mapped effect-allele column is a silent sign flip, not an error, so this runs
first and the config is corrected to match it.


::: {.cell}

```{.r .cell-code}
sources <- c(
  purrr::imap(cfg$exposures, ~ list(label = .x$label, file = .x$file)),
  purrr::imap(cfg$outcomes,  ~ list(label = .x$label, file = .x$file))
)

for (s in sources) {
  p <- pp(s$file)
  cat("\n----", s$label, "----\n")
  if (!file.exists(p)) {
    cat("  MISSING:", s$file, "\n  Run: bash scripts/fetch_robust_mr_data.sh\n")
  } else {
    cat(inspect_gwas_header(p, n = 2), sep = "\n")
  }
}
```

::: {.cell-output .cell-output-stdout}

```

---- LDL-C (GLGC 2021, EUR) ----
rsID	CHROM	POS_b37	REF	ALT	N	N_studies	POOLED_ALT_AF	EFFECT_SIZE	SE	pvalue_neg_log10	pvalue	pvalue_neg_log10_GC	pvalue_GC
rs367896724	1	    10177	A	AC	7977	6	0.349	-0.0369442	0.027669	0.740395661246388	0.182	0.603293115674572	0.249

---- LDL-C (GLGC 2013, pre-UKB) ----
SNP_hg18	SNP_hg19	rsid	A1	A2	beta	se	N	P-value	Freq.A1.1000G.EUR
chr10:10000135	chr10:9960129	rs4747841	a	g	0.0037	0.0052	89138.00	0.7158	0.4908

---- Total cholesterol (GLGC 2021, EUR) ----
rsID	CHROM	POS_b37	REF	ALT	N	N_studies	POOLED_ALT_AF	EFFECT_SIZE	SE	pvalue_neg_log10	pvalue	pvalue_neg_log10_GC	pvalue_GC
rs367896724	1	    10177	A	AC	8577	5	0.349	-0.0146697	0.0266621	0.234945153469547	0.582	0.189364074110782	0.647

---- Heel eBMD (Morris 2019, UKB) ----
variant_id	snp.1	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	info	beta	standard_error	p_value	n	ci_upper	odds_ratio	ci_lower
1:55326:T:C	rs3107975	1	55326	T	C	0.99165	0.313908	-0.0324378	0.0181701	2.5E-01	1.3E-01	NA	NA	NA
```


:::
:::


## Step 2 — read, QC, and harmonise

QC follows the MR-APSS published criteria so that both methods see an identical
variant set: HapMap3 restriction, MAF > 0.05, INFO > 0.9, unambiguous biallelic
alleles, MHC excluded, and a $\chi^2$ cap of max(N/1000, 80).


::: {.cell}

```{.r .cell-code}
hm3 <- read_hm3(pp(cfg$paths$hm3))

prep_one <- function(spec, label) {
  p <- pp(spec$file)
  if (!file.exists(p)) return(NULL)
  raw <- if (isTRUE(spec$vcf)) {
    read_gwasvcf(p, n_fixed = spec$n_fixed)          # OpenGWAS bulk GWAS-VCF
  } else {
    read_gwas(p, cols = spec$cols, n_fixed = spec$n_fixed, log10p = spec$log10p)
  }
  qc_gwas(raw, hm3 = hm3,
          maf_min  = cfg$qc$maf_min,
          info_min = cfg$qc$info_min,
          drop_mhc = cfg$qc$drop_mhc,
          label    = label)
}

gwas <- list(
  glgc2021_ldl   = prep_one(cfg$exposures$glgc2021_ldl,   "LDL (GLGC 2021)"),
  willer2013_ldl = prep_one(cfg$exposures$willer2013_ldl, "LDL (GLGC 2013)"),
  ebmd           = prep_one(cfg$outcomes$ebmd_morris2019, "eBMD (Morris 2019)")
) |> purrr::compact()

qc_table <- purrr::map_dfr(gwas, qc_log)
qc_table |>
  pivot_wider(names_from = trait, values_from = n) |>
  kable(caption = "QC attrition by step", format.args = list(big.mark = ",")) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>QC attrition by step</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> step </th>
   <th style="text-align:right;"> LDL (GLGC 2021) </th>
   <th style="text-align:right;"> LDL (GLGC 2013) </th>
   <th style="text-align:right;"> eBMD (Morris 2019) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> input </td>
   <td style="text-align:right;"> 47,006,483 </td>
   <td style="text-align:right;"> 2,437,751 </td>
   <td style="text-align:right;"> 13,753,401 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> unique rsID </td>
   <td style="text-align:right;"> 44,473,875 </td>
   <td style="text-align:right;"> 2,436,957 </td>
   <td style="text-align:right;"> 13,753,401 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> biallelic ACGT </td>
   <td style="text-align:right;"> 41,556,903 </td>
   <td style="text-align:right;"> 2,436,957 </td>
   <td style="text-align:right;"> 13,753,343 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> unambiguous strand </td>
   <td style="text-align:right;"> 35,282,762 </td>
   <td style="text-align:right;"> 2,061,312 </td>
   <td style="text-align:right;"> 11,692,171 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> in HapMap3 </td>
   <td style="text-align:right;"> 1,215,630 </td>
   <td style="text-align:right;"> 1,042,744 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAF &gt; 0.05 </td>
   <td style="text-align:right;"> 1,073,848 </td>
   <td style="text-align:right;"> 960,847 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MHC removed </td>
   <td style="text-align:right;"> 1,071,860 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> chi2 &lt; 1320 </td>
   <td style="text-align:right;"> 1,071,731 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> chi2 &lt; 90 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 960,452 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> chi2 &lt; NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 3 — write the two input formats

CAUSE and MR-APSS want different shapes. Both are written from the same QC'd
object so they cannot silently diverge.


::: {.cell}

```{.r .cell-code}
dir.create(pp(cfg$paths$cache), recursive = TRUE, showWarnings = FALSE)

walk2(gwas, names(gwas), function(d, nm) {
  saveRDS(to_apss(d),  pp(cfg$paths$cache, paste0(nm, "_apss.rds")))
  saveRDS(to_cause(d), pp(cfg$paths$cache, paste0(nm, "_cause.rds")))
})

tibble(
  dataset = names(gwas),
  n_snps  = map_int(gwas, nrow),
  median_N = map_dbl(gwas, ~ median(.x$n, na.rm = TRUE)),
  mean_chi2 = map_dbl(gwas, ~ mean(.x$z^2, na.rm = TRUE))
) |>
  kable(caption = "Analysis-ready datasets. Mean chi-square is the LDSC power diagnostic; below ~1.02 the LDSC-based background model will be unstable.",
        digits = 3, format.args = list(big.mark = ",")) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Analysis-ready datasets. Mean chi-square is the LDSC power diagnostic; below ~1.02 the LDSC-based background model will be unstable.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> dataset </th>
   <th style="text-align:right;"> n_snps </th>
   <th style="text-align:right;"> median_N </th>
   <th style="text-align:right;"> mean_chi2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> glgc2021_ldl </td>
   <td style="text-align:right;"> 1,071,731 </td>
   <td style="text-align:right;"> 1,320,016 </td>
   <td style="text-align:right;"> 3.187 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> willer2013_ldl </td>
   <td style="text-align:right;"> 960,452 </td>
   <td style="text-align:right;"> 89,872 </td>
   <td style="text-align:right;"> 1.217 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ebmd </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NaN </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 4 — sanity check against the existing estimates

Before running anything Bayesian, confirm that a plain IVW on this variant set
reproduces the effect already reported in `ANALYSIS.md` (all-cholesterol → heel
BMD, β ≈ −0.051). If it does not, the harmonisation is wrong and nothing
downstream is interpretable.


::: {.cell}

```{.r .cell-code}
ivw_check <- function(exp_name) {
  e <- gwas[[exp_name]]; o <- gwas$ebmd
  if (is.null(e) || is.null(o)) return(NULL)

  iv <- e |>
    filter(pval < 5e-8) |>
    clump_local(snp_col = "snp", p_col = "pval",
                r2 = 0.001, kb = 1000, p_thresh = 5e-8,
                bfile = pp(cfg$paths$plink_bfile),
                plink_bin = cfg$paths$plink_bin)

  d <- iv |>
    select(snp, a1, a2, b_exp = beta, se_exp = se) |>
    inner_join(o |> select(snp, a1_o = a1, a2_o = a2, b_out = beta, se_out = se),
               by = "snp") |>
    mutate(b_out = if_else(a1 == a1_o, b_out,
                           if_else(a1 == a2_o, -b_out, NA_real_))) |>
    filter(!is.na(b_out))

  w <- 1 / d$se_out^2
  b <- sum(w * d$b_exp * d$b_out) / sum(w * d$b_exp^2)
  se <- sqrt(1 / sum(w * d$b_exp^2))
  tibble(exposure = exp_name, n_iv = nrow(d), b = b, se = se,
         pval = 2 * pnorm(abs(b / se), lower.tail = FALSE),
         f_mean = mean((d$b_exp / d$se_exp)^2))
}

ivw_baseline <- map_dfr(c("glgc2021_ldl", "willer2013_ldl"), ivw_check)
write_csv(ivw_baseline, pp(cfg$paths$out, "ivw_baseline.csv"))
ivw_baseline |>
  kable(caption = "Conventional IVW on the same harmonised data (sanity anchor). Mean F should be well above 10; compare b to the -0.051 in ANALYSIS.md.",
        digits = 4) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Conventional IVW on the same harmonised data (sanity anchor). Mean F should be well above 10; compare b to the -0.051 in ANALYSIS.md.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> exposure </th>
   <th style="text-align:right;"> n_iv </th>
   <th style="text-align:right;"> b </th>
   <th style="text-align:right;"> se </th>
   <th style="text-align:right;"> pval </th>
   <th style="text-align:right;"> f_mean </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> glgc2021_ldl </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> NaN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> willer2013_ldl </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> NaN </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Next

- `robust_mr_apss.qmd` — MR-APSS. Fast (minutes); run locally.
- `robust_mr_cause.qmd` — CAUSE. Memory-hungry; see the Great Lakes script.
- `robust_mr_summary.qmd` — reconciles everything against the conventional estimates.
