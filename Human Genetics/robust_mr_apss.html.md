---
title: "Robust MR — MR-APSS (pleiotropy + sample structure + winner's curse)"
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

MR-APSS is the single method in this project that addresses all three of our
concerns simultaneously. It decomposes each observed SNP effect into a
**background** component — polygenicity ($\Omega$), correlated pleiotropy, and
sample structure ($C$) — and a **foreground** component, on which causal
inference is performed while allowing uncorrelated pleiotropy.

Three properties matter for us:

1. **$C$ absorbs sample overlap.** $C$ is a 2×2 matrix of LDSC intercepts; the
   diagonal captures inflation within each GWAS, and $C_{12}$ — the cross-trait
   intercept — is proportional to the sample overlap times the phenotypic
   correlation. This makes overlap a *measured quantity*, not an assumption.
2. **Relaxed instrument threshold with an explicit selection-bias correction.**
   IVs are taken at $p<5\times10^{-5}$ rather than $5\times10^{-8}$. On its own
   that would import winner's curse and weak-instrument bias; `Cor.SelectionBias
   = TRUE` corrects for the selection step, which is what makes the relaxed
   threshold a power gain rather than a bias.
3. **$\Omega$ soaks up correlated pleiotropy** before the causal parameter is
   estimated, so a shared polygenic background between lipids and bone does not
   masquerade as a causal effect.

Runtime is minutes, not hours — this one does not need Great Lakes.


::: {.cell cache.extra='["300fa520e3962a0596bb81c68cead500","409821511791fbd067ed06f24f5faeef"]'}

```{.r .cell-code}
read_cached <- function(nm) {
  p <- pp(cfg$paths$cache, paste0(nm, "_apss.rds"))
  if (!file.exists(p)) stop("Run robust_mr_prep.qmd first — missing ", p)
  readRDS(p)
}

dat <- list(
  glgc2021_ldl   = read_cached("glgc2021_ldl"),
  willer2013_ldl = read_cached("willer2013_ldl"),
  ebmd           = read_cached("ebmd")
)

map_dfr(dat, ~ tibble(n_snp = nrow(.x), median_N = median(.x$N)), .id = "dataset") |>
  kable(caption = "MR-APSS input datasets (SNP / A1 / A2 / Z / P / N)",
        format.args = list(big.mark = ",")) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>MR-APSS input datasets (SNP / A1 / A2 / Z / P / N)</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> dataset </th>
   <th style="text-align:right;"> n_snp </th>
   <th style="text-align:right;"> median_N </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> glgc2021_ldl </td>
   <td style="text-align:right;"> 1,071,731 </td>
   <td style="text-align:right;"> 1,320,016 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> willer2013_ldl </td>
   <td style="text-align:right;"> 960,452 </td>
   <td style="text-align:right;"> 89,872 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ebmd </td>
   <td style="text-align:right;"> 1,030,710 </td>
   <td style="text-align:right;"> 426,824 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 1 — background parameters ($C$, $\Omega$) by bivariate LDSC


::: {.cell cache.extra='["300fa520e3962a0596bb81c68cead500","409821511791fbd067ed06f24f5faeef"]'}

```{.r .cell-code}
library(MRAPSS)

ldsc_dir <- pp(cfg$paths$ldscore)
stopifnot("LD score directory not found — see scripts/fetch_robust_mr_data.sh" =
            dir.exists(ldsc_dir))

# Read once, pass as ld/M rather than ldscore.dir: est_paras() would otherwise
# re-read all 22 chromosomes per trait pair, and its hard-coded filename pattern
# only matches the eur_w_ld_chr layout. read_ldscores() handles both layouts.
ldsc <- read_ldscores(ldsc_dir)
```

::: {.cell-output .cell-output-stderr}

```
Read LD scores for 1,190,321 SNPs (M = 5,961,159)
```


:::

```{.r .cell-code}
fit_paras <- function(exp_name) {
  MRAPSS::est_paras(
    dat1 = dat[[exp_name]],
    dat2 = dat$ebmd,
    trait1.name = cfg$exposures[[exp_name]]$label,
    trait2.name = cfg$outcomes$ebmd_morris2019$label,
    ld = ldsc$ld,
    M  = ldsc$M
  )
}

paras <- list(
  overlapping  = fit_paras("glgc2021_ldl"),
  `overlap-free` = fit_paras("willer2013_ldl")
)
```

::: {.cell-output .cell-output-stderr}

```
Merge dat1 and dat2 by SNP ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Harmonize the direction of SNP effects of exposure and outcome
```


:::

::: {.cell-output .cell-output-stderr}

```
Read in LD scores ... 
```


:::

::: {.cell-output .cell-output-stderr}

```
Add LD scores to the harmonized data set...
```


:::

::: {.cell-output .cell-output-stderr}

```
The Harmonized data set will also be used for MR analysis 
```


:::

::: {.cell-output .cell-output-stderr}

```
Begin estimation of C and Omega using LDSC ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Estimate heritability for trait 1 ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Using two-step estimator with cutoff at 30.
```


:::

::: {.cell-output .cell-output-stderr}

```
Mean Chi2:3.2221.
```


:::

::: {.cell-output .cell-output-stderr}

```
Intercept: 1.1604(0.0252).
```


:::

::: {.cell-output .cell-output-stderr}

```
Total Observed scale h2:0.0845(0.0096).
```


:::

::: {.cell-output .cell-output-stderr}

```
Estimate heritability for trait 2 ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Using two-step estimator with cutoff at 30.
```


:::

::: {.cell-output .cell-output-stderr}

```
Mean Chi2:4.3766.
```


:::

::: {.cell-output .cell-output-stderr}

```
Intercept: 1.5137 (0.0374).
```


:::

::: {.cell-output .cell-output-stderr}

```
Total Observed scale h2:0.3395 (0.0185).
```


:::

::: {.cell-output .cell-output-stderr}

```
Estimate genetic covariance ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Using two-step estimator with cutoff at 30.
```


:::

::: {.cell-output .cell-output-stderr}

```
Intercept: 0.0226 (0.0122).
```


:::

::: {.cell-output .cell-output-stderr}

```
Total Observed scale gencov: -0.0063 (0.0026).
```


:::

::: {.cell-output .cell-output-stderr}

```
Merge dat1 and dat2 by SNP ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Harmonize the direction of SNP effects of exposure and outcome
```


:::

::: {.cell-output .cell-output-stderr}

```
Read in LD scores ... 
```


:::

::: {.cell-output .cell-output-stderr}

```
Add LD scores to the harmonized data set...
```


:::

::: {.cell-output .cell-output-stderr}

```
The Harmonized data set will also be used for MR analysis 
```


:::

::: {.cell-output .cell-output-stderr}

```
Begin estimation of C and Omega using LDSC ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Estimate heritability for trait 1 ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Using two-step estimator with cutoff at 30.
```


:::

::: {.cell-output .cell-output-stderr}

```
Mean Chi2:1.2186.
```


:::

::: {.cell-output .cell-output-stderr}

```
Intercept: 1.0011(0.0108).
```


:::

::: {.cell-output .cell-output-stderr}

```
Total Observed scale h2:0.1192(0.012).
```


:::

::: {.cell-output .cell-output-stderr}

```
Estimate heritability for trait 2 ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Using two-step estimator with cutoff at 30.
```


:::

::: {.cell-output .cell-output-stderr}

```
Mean Chi2:4.3339.
```


:::

::: {.cell-output .cell-output-stderr}

```
Intercept: 1.5176 (0.0383).
```


:::

::: {.cell-output .cell-output-stderr}

```
Total Observed scale h2:0.3343 (0.0188).
```


:::

::: {.cell-output .cell-output-stderr}

```
Estimate genetic covariance ...
```


:::

::: {.cell-output .cell-output-stderr}

```
Using two-step estimator with cutoff at 30.
```


:::

::: {.cell-output .cell-output-stderr}

```
Intercept: -0.0074 (0.0095).
```


:::

::: {.cell-output .cell-output-stderr}

```
Total Observed scale gencov: -0.0056 (0.0046).
```


:::

```{.r .cell-code}
saveRDS(paras, pp(cfg$paths$cache, "apss_paras.rds"))
```
:::


### The overlap measurement

This is the table that answers the question directly. $C_{12}$ should be
materially non-zero for the GLGC 2021 arm (UK Biobank on both sides) and close to
zero for the pre-UK-Biobank GLGC 2013 arm. If it is *not* near zero in the second
arm, either the arm is not as clean as we think or there is shared population
stratification rather than shared samples.


::: {.cell}

```{.r .cell-code}
bg <- imap_dfr(paras, ~ tidy_background(
  .x,
  exposure = cfg$exposures[[if (.y == "overlapping") "glgc2021_ldl" else "willer2013_ldl"]]$label,
  outcome  = cfg$outcomes$ebmd_morris2019$label,
  arm      = .y
))
write_csv(bg, pp(cfg$paths$out, "apss_background_parameters.csv"))

bg |>
  select(arm, C11, C22, C12, overlap_flag, rg) |>
  kable(caption = "MR-APSS background parameters. C11/C22 are the univariate LDSC intercepts (confounding within each GWAS); C12 is the cross-trait intercept and is the direct estimate of sample overlap. rg is the genetic correlation implied by Omega — the correlated-pleiotropy channel.",
        digits = 4) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>MR-APSS background parameters. C11/C22 are the univariate LDSC intercepts (confounding within each GWAS); C12 is the cross-trait intercept and is the direct estimate of sample overlap. rg is the genetic correlation implied by Omega — the correlated-pleiotropy channel.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> C11 </th>
   <th style="text-align:right;"> C22 </th>
   <th style="text-align:right;"> C12 </th>
   <th style="text-align:left;"> overlap_flag </th>
   <th style="text-align:right;"> rg </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> 1.1604 </td>
   <td style="text-align:right;"> 1.5137 </td>
   <td style="text-align:right;"> 0.0226 </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:right;"> -0.0371 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> 1.0011 </td>
   <td style="text-align:right;"> 1.5176 </td>
   <td style="text-align:right;"> -0.0074 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> -0.0280 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 2 — IV selection and LD clumping


::: {.cell cache.extra='["300fa520e3962a0596bb81c68cead500","409821511791fbd067ed06f24f5faeef"]'}

```{.r .cell-code}
bfile     <- pp(cfg$paths$plink_bfile)
plink_bin <- cfg$paths$plink_bin

# MRAPSS ships its own clump() wrapper, but it has moved between versions.
# Fall back to the project's ieugwasr-based clumper, which produces the same thing.
clump_arm <- function(par) {
  if (!is.null(getNamespace("MRAPSS")$clump)) {
    out <- try(MRAPSS::clump(
      par$dat,
      IV.Threshold = cfg$mrapss$iv_threshold,
      SNP_col   = "SNP",
      pval_col  = "pval.exp",
      clump_kb  = cfg$mrapss$clump_kb,
      clump_r2  = cfg$mrapss$clump_r2,
      bfile     = bfile,
      plink_bin = plink_bin
    ), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    message("MRAPSS::clump() failed; falling back to clump_local().")
  }
  clump_local(par$dat, snp_col = "SNP", p_col = "pval.exp",
              r2 = cfg$mrapss$clump_r2, kb = cfg$mrapss$clump_kb,
              p_thresh = cfg$mrapss$iv_threshold,
              bfile = bfile, plink_bin = plink_bin)
}

MRdat <- map(paras, clump_arm) |>
  map(sanitise_mrdat, iv_threshold = cfg$mrapss$iv_threshold)
saveRDS(MRdat, pp(cfg$paths$cache, "apss_MRdat.rds"))

map_dfr(MRdat, ~ tibble(
  n_iv         = nrow(.x),
  threshold    = unique(.x$Threshold),
  min_pval_exp = min(.x$pval.exp),
  max_pval_exp = max(.x$pval.exp)
), .id = "arm") |>
  kable(caption = sprintf("Instruments after selection at p < %.0e and clumping at r2 = %.3f. max_pval_exp must not exceed the threshold — if it does, the clumping step did not actually filter on p.",
                          cfg$mrapss$iv_threshold, cfg$mrapss$clump_r2)) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Instruments after selection at p &lt; 5e-05 and clumping at r2 = 0.001. max_pval_exp must not exceed the threshold — if it does, the clumping step did not actually filter on p.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> n_iv </th>
   <th style="text-align:right;"> threshold </th>
   <th style="text-align:right;"> min_pval_exp </th>
   <th style="text-align:right;"> max_pval_exp </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> 948 </td>
   <td style="text-align:right;"> 5e-05 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4.90e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> 164 </td>
   <td style="text-align:right;"> 5e-05 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4.85e-05 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


::: callout-note
The instrument count here should be in the hundreds to low thousands — far more
than the ~230–280 SNPs the genome-wide-significant analyses in this repo use. That
is the point of the 5e-5 threshold, and it is only defensible because
`Cor.SelectionBias = TRUE` corrects the resulting winner's curse.
:::

## Step 3 — fit MR-APSS

Each arm is fitted twice: once with the estimated $C$ (the correction switched on)
and once with $C = I$ (switched off). The gap between them *is* the sample-overlap
bias, expressed in the units of the outcome.


::: {.cell cache.extra='["300fa520e3962a0596bb81c68cead500","409821511791fbd067ed06f24f5faeef"]'}

```{.r .cell-code}
exp_of <- function(arm) if (arm == "overlapping") "glgc2021_ldl" else "willer2013_ldl"

fit_apss <- function(arm, use_C) {
  par <- paras[[arm]]
  MRAPSS::MRAPSS(
    MRdat[[arm]],
    exposure = cfg$exposures[[exp_of(arm)]]$label,
    outcome  = cfg$outcomes$ebmd_morris2019$label,
    C        = if (use_C) par$C else diag(2),
    Omega    = par$Omega,
    Cor.SelectionBias = cfg$mrapss$cor_selection_bias
  )
}

grid <- tidyr::expand_grid(arm = names(paras), use_C = c(TRUE, FALSE))
if (!isTRUE(cfg$mrapss$run_C_identity)) grid <- filter(grid, use_C)

fits <- pmap(grid, fit_apss)
```

::: {.cell-output .cell-output-stdout}

```
***********************************************************
MR test results of  LDL-C (GLGC 2021, EUR)  on  Heel eBMD (Morris 2019, UKB) : 
MR-APSS: beta =  -0.0561 , beta.se =  0.0265 , p-value =  3.4388e-02 . 
Total NO. of IVs=  948 , NO. of valid IVs with foreground signals:  412.0344 . 
***********************************************************
***********************************************************
MR test results of  LDL-C (GLGC 2021, EUR)  on  Heel eBMD (Morris 2019, UKB) : 
MR-APSS: beta =  -0.0517 , beta.se =  0.0262 , p-value =  4.8075e-02 . 
Total NO. of IVs=  948 , NO. of valid IVs with foreground signals:  427.4942 . 
***********************************************************
***********************************************************
MR test results of  LDL-C (GLGC 2013, pre-UKB)  on  Heel eBMD (Morris 2019, UKB) : 
MR-APSS: beta =  -0.0375 , beta.se =  0.0187 , p-value =  4.4928e-02 . 
Total NO. of IVs=  164 , NO. of valid IVs with foreground signals:  129.2338 . 
***********************************************************
***********************************************************
MR test results of  LDL-C (GLGC 2013, pre-UKB)  on  Heel eBMD (Morris 2019, UKB) : 
MR-APSS: beta =  -0.0361 , beta.se =  0.0184 , p-value =  4.9469e-02 . 
Total NO. of IVs=  164 , NO. of valid IVs with foreground signals:  129.8323 . 
***********************************************************
```


:::

```{.r .cell-code}
saveRDS(list(grid = grid, fits = fits), pp(cfg$paths$cache, "apss_fits.rds"))

apss_res <- pmap_dfr(list(grid$arm, grid$use_C, fits), function(arm, use_C, fit) {
  tidy_apss(fit,
            exposure  = cfg$exposures[[exp_of(arm)]]$label,
            outcome   = cfg$outcomes$ebmd_morris2019$label,
            arm       = arm,
            C_setting = if (use_C) "corrected" else "uncorrected")
})
write_csv(apss_res, pp(cfg$paths$out, "mrapss_results.csv"))

apss_res |>
  select(arm, C_setting, b, se, lci, uci, pval, n_iv, n_valid_iv, pi0) |>
  kable(caption = "MR-APSS causal effect of LDL-C on heel eBMD (SD per SD). 'uncorrected' rows are the identical model refitted with C = I, i.e. the sample-structure correction disabled. n_valid_iv is the posterior-expected number of instruments carrying foreground signal, sum(post$Pi).",
        digits = 4) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>MR-APSS causal effect of LDL-C on heel eBMD (SD per SD). 'uncorrected' rows are the identical model refitted with C = I, i.e. the sample-structure correction disabled. n_valid_iv is the posterior-expected number of instruments carrying foreground signal, sum(post$Pi).</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:left;"> C_setting </th>
   <th style="text-align:right;"> b </th>
   <th style="text-align:right;"> se </th>
   <th style="text-align:right;"> lci </th>
   <th style="text-align:right;"> uci </th>
   <th style="text-align:right;"> pval </th>
   <th style="text-align:right;"> n_iv </th>
   <th style="text-align:right;"> n_valid_iv </th>
   <th style="text-align:right;"> pi0 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> corrected </td>
   <td style="text-align:right;"> -0.0561 </td>
   <td style="text-align:right;"> 0.0265 </td>
   <td style="text-align:right;"> -0.1080 </td>
   <td style="text-align:right;"> -0.0041 </td>
   <td style="text-align:right;"> 0.0344 </td>
   <td style="text-align:right;"> 948 </td>
   <td style="text-align:right;"> 412.0364 </td>
   <td style="text-align:right;"> 0.4346 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> uncorrected </td>
   <td style="text-align:right;"> -0.0517 </td>
   <td style="text-align:right;"> 0.0262 </td>
   <td style="text-align:right;"> -0.1030 </td>
   <td style="text-align:right;"> -0.0004 </td>
   <td style="text-align:right;"> 0.0481 </td>
   <td style="text-align:right;"> 948 </td>
   <td style="text-align:right;"> 427.4945 </td>
   <td style="text-align:right;"> 0.4509 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> corrected </td>
   <td style="text-align:right;"> -0.0375 </td>
   <td style="text-align:right;"> 0.0187 </td>
   <td style="text-align:right;"> -0.0741 </td>
   <td style="text-align:right;"> -0.0008 </td>
   <td style="text-align:right;"> 0.0449 </td>
   <td style="text-align:right;"> 164 </td>
   <td style="text-align:right;"> 129.2384 </td>
   <td style="text-align:right;"> 0.7880 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> uncorrected </td>
   <td style="text-align:right;"> -0.0361 </td>
   <td style="text-align:right;"> 0.0184 </td>
   <td style="text-align:right;"> -0.0721 </td>
   <td style="text-align:right;"> -0.0001 </td>
   <td style="text-align:right;"> 0.0495 </td>
   <td style="text-align:right;"> 164 </td>
   <td style="text-align:right;"> 129.8373 </td>
   <td style="text-align:right;"> 0.7917 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


### How much did the overlap correction move the estimate?


::: {.cell}

```{.r .cell-code}
if (n_distinct(apss_res$C_setting) == 2) {
  delta <- apss_res |>
    select(arm, C_setting, b, se) |>
    pivot_wider(names_from = C_setting, values_from = c(b, se)) |>
    mutate(shift       = b_corrected - b_uncorrected,
           shift_in_se = shift / se_corrected)
  write_csv(delta, pp(cfg$paths$out, "mrapss_overlap_shift.csv"))
  delta |>
    kable(caption = "Effect of switching the sample-structure correction on. A shift of more than ~0.5 SE in the overlapping arm but not the overlap-free arm is the signature of genuine overlap bias.",
          digits = 4) |>
    kable_styling(full_width = FALSE)
}
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Effect of switching the sample-structure correction on. A shift of more than ~0.5 SE in the overlapping arm but not the overlap-free arm is the signature of genuine overlap bias.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> b_corrected </th>
   <th style="text-align:right;"> b_uncorrected </th>
   <th style="text-align:right;"> se_corrected </th>
   <th style="text-align:right;"> se_uncorrected </th>
   <th style="text-align:right;"> shift </th>
   <th style="text-align:right;"> shift_in_se </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0561 </td>
   <td style="text-align:right;"> -0.0517 </td>
   <td style="text-align:right;"> 0.0265 </td>
   <td style="text-align:right;"> 0.0262 </td>
   <td style="text-align:right;"> -0.0044 </td>
   <td style="text-align:right;"> -0.1655 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0375 </td>
   <td style="text-align:right;"> -0.0361 </td>
   <td style="text-align:right;"> 0.0187 </td>
   <td style="text-align:right;"> 0.0184 </td>
   <td style="text-align:right;"> -0.0014 </td>
   <td style="text-align:right;"> -0.0749 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 4 — diagnostic plot


::: {.cell}

```{.r .cell-code}
walk(seq_along(fits), function(i) {
  ok <- try(print(MRAPSS::MRplot(
    fits[[i]],
    exposure = "LDL-C",
    outcome  = "heel eBMD"
  )), silent = TRUE)
  invisible(ok)
})
```

::: {.cell-output-display}
![](robust_mr_apss_files/figure-html/plot-1.png){width=576}
:::

::: {.cell-output-display}
![](robust_mr_apss_files/figure-html/plot-2.png){width=576}
:::

::: {.cell-output-display}
![](robust_mr_apss_files/figure-html/plot-3.png){width=576}
:::

::: {.cell-output-display}
![](robust_mr_apss_files/figure-html/plot-4.png){width=576}
:::
:::


Triangles are observed SNP effects; colour is the posterior probability that the
instrument carries foreground (i.e. valid) signal — the per-SNP `post$Pi`, whose
sum is the `n_valid_iv` column above. Expect it to be well below the raw
instrument count: that gap is the model discounting background-driven variants,
and it is the mechanism by which the relaxed $5\times10^{-5}$ threshold buys power
without buying bias.

## Interpretation notes

- The MR-APSS $\beta$ is on a **per-SD-exposure, per-SD-outcome** scale by
  construction ($b = Z/\sqrt{N}$ internally), so it is directly comparable to the
  IVW estimates in `ANALYSIS.md` provided those were also run on standardised
  traits.
- A large $C_{11}$ or $C_{22}$ (LDSC intercept well above 1) means one of the
  GWAS is inflated by confounding. The MR-APSS authors flag that in that regime
  the causal estimate carries larger standard errors than IVW or RAPS — that is
  the correction being honest, not a failure.
- MR-APSS assumes LDSC's assumptions hold, which means it inherits LDSC's
  requirements: European-ancestry samples matched to `eur_w_ld_chr`, N > ~5k, and
  HapMap3 variants. All satisfied here.
