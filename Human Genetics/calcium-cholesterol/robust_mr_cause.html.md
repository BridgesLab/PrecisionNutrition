---
title: "Robust MR — CAUSE (correlated horizontal pleiotropy)"
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

CAUSE tests whether the LDL → eBMD signal is better explained by a **causal
effect** or by a **shared factor** $U$ acting on both traits — correlated
horizontal pleiotropy. This is the failure mode that Egger, weighted median and
mode-based estimators cannot detect, because a shared factor produces exactly the
proportional, low-heterogeneity pattern those methods read as a clean causal
effect.

Two structural differences from every other MR in this repo:

1. **It is a model comparison, not an effect estimate.** The primary output is an
   ELPD difference between the sharing and causal models. $\hat\gamma$ exists, but
   the inferential claim lives in the comparison. CAUSE is deliberately
   conservative; a null here means "a shared factor explains this as well as
   causation does", not "no effect".
2. **It is robust to sample overlap** because the nuisance parameter $\rho$ — the
   correlation of test statistics under the null — is estimated genome-wide from
   ~1M variants, and overlap is one of the things that induces exactly that
   correlation.

::: callout-important
## Runtime and memory

`est_cause_params()` plus the three-model fit is the expensive part: expect
**1–2 hours and 16–32 GB** per arm on the HapMap3-restricted input built by
`robust_mr_prep.qmd` (~1M variants — which is also exactly the number of variants
CAUSE wants for parameter estimation, so nothing is lost by the restriction).

That is laptop-feasible overnight, but `scripts/run_cause_greatlakes.R` and
`scripts/cause_greatlakes.sbatch` run the chunks below headless so both arms can go
in parallel on Great Lakes. Either way the result lands in
`cause_fit_<arm>.rds`, and the chunk below picks it up and skips recomputation.
:::

## Step 1 — merge and estimate nuisance parameters


::: {.cell cache.extra='["300fa520e3962a0596bb81c68cead500","b1b592161a47e5212d1889df403e39c7"]'}

```{.r .cell-code}
library(cause)

cache_rds <- function(nm) pp(cfg$paths$cache, nm)

run_or_load <- function(arm, exp_name) {
  fit_path <- cache_rds(paste0("cause_fit_", arm, ".rds"))
  if (file.exists(fit_path)) {
    message("Loading precomputed CAUSE fit for arm: ", arm)
    return(readRDS(fit_path))
  }

  message("No cached fit for '", arm, "'. Computing locally — this is slow. ",
          "Consider: sbatch scripts/cause_greatlakes.sbatch ", arm)

  e <- readRDS(cache_rds(paste0(exp_name, "_cause.rds")))
  o <- readRDS(cache_rds("ebmd_cause.rds"))

  X <- cause::gwas_merge(
    e, o,
    snp_name_cols = c("SNP", "SNP"),
    beta_hat_cols = c("beta_hat", "beta_hat"),
    se_cols       = c("se", "se"),
    A1_cols       = c("A1", "A1"),
    A2_cols       = c("A2", "A2")
  )

  varlist <- with(X, sample(snp,
                            size = min(cfg$cause$n_param_variants, nrow(X)),
                            replace = FALSE))
  params <- cause::est_cause_params(X, varlist)

  # LD pruning: p < 1e-3, r2 < 0.01, using the project's local plink panel.
  # (cause::ld_prune() with the Zenodo 1000G blocks is the alternative; see below.)
  X <- X |> mutate(pval1 = 2 * pnorm(abs(beta_hat_1 / seb1), lower.tail = FALSE))
  top <- clump_local(X, snp_col = "snp", p_col = "pval1",
                     r2 = cfg$cause$prune_r2, kb = cfg$cause$prune_kb,
                     p_thresh = cfg$cause$prune_p,
                     bfile = pp(cfg$paths$plink_bfile),
                     plink_bin = cfg$paths$plink_bin)

  fit <- cause::cause(X = X, variants = top$snp, param_ests = params)
  out <- list(fit = fit, params = params, n_variants = length(top$snp), n_merged = nrow(X))
  saveRDS(out, fit_path)
  out
}

cause_out <- list(
  overlapping    = run_or_load("overlapping",   "glgc2021_ldl"),
  `overlap-free` = run_or_load("overlap-free",  "willer2013_ldl")
)
```

::: {.cell-output .cell-output-stderr}

```
Loading precomputed CAUSE fit for arm: overlapping
```


:::

::: {.cell-output .cell-output-stderr}

```
Loading precomputed CAUSE fit for arm: overlap-free
```


:::
:::


::: callout-note
## Alternative pruning with CAUSE's own LD blocks

If you would rather use the reference LD the CAUSE authors ship (1000 Genomes EUR,
Zenodo record 1464357) instead of plink clumping:

```r
pruned <- purrr::map_dfr(1:22, function(chr) {
  ld  <- readRDS(file.path(cfg$paths$cause_ld, sprintf("chr%s_AF0.05_0.1.RDS", chr)))
  snp <- readRDS(file.path(cfg$paths$cause_ld, sprintf("chr%s_AF0.05_snpdata.RDS", chr)))
  cause::ld_prune(variants = X |> dplyr::filter(snp %in% snp$SNP),
                  ld = ld, total_ld_variants = snp$SNP,
                  pval_cols = c("pval1"), pval_thresh = c(cfg$cause$prune_p))
})
```

The two give very similar instrument sets; plink clumping is used by default
because the panel is already in this repo for the MVMR work.
:::

### Nuisance parameters

$\rho$ is the correlation of the two traits' test statistics under the null. It is
non-zero when the GWAS share samples, share stratification, or the traits are
genuinely correlated for non-causal reasons — and it is the reason CAUSE tolerates
overlap.


::: {.cell}

```{.r .cell-code}
rho_tab <- imap_dfr(cause_out, ~ tibble(
  arm = .y,
  rho = .x$params$rho,
  n_merged   = .x$n_merged,
  n_variants = .x$n_variants
))
write_csv(rho_tab, pp(cfg$paths$out, "cause_nuisance.csv"))
rho_tab |>
  kable(caption = "CAUSE nuisance parameters. Larger |rho| in the overlapping arm than the overlap-free arm is the expected signature of shared samples.",
        digits = 4, format.args = list(big.mark = ",")) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>CAUSE nuisance parameters. Larger |rho| in the overlapping arm than the overlap-free arm is the expected signature of shared samples.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> rho </th>
   <th style="text-align:right;"> n_merged </th>
   <th style="text-align:right;"> n_variants </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0109 </td>
   <td style="text-align:right;"> 1,025,362 </td>
   <td style="text-align:right;"> 2,089 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0146 </td>
   <td style="text-align:right;"> 909,826 </td>
   <td style="text-align:right;"> 651 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 2 — model comparison


::: {.cell}

```{.r .cell-code}
elpd_tab <- imap_dfr(cause_out, function(x, arm) {
  as_tibble(x$fit$elpd) |> mutate(arm = arm, .before = 1)
})
write_csv(elpd_tab, pp(cfg$paths$out, "cause_elpd.csv"))
elpd_tab |>
  kable(caption = "ELPD model comparison. The decisive row is model1 = sharing vs model2 = causal: a negative delta_elpd with a significant z favours the causal model.",
        digits = 3) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>ELPD model comparison. The decisive row is model1 = sharing vs model2 = causal: a negative delta_elpd with a significant z favours the causal model.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:left;"> model1 </th>
   <th style="text-align:left;"> model2 </th>
   <th style="text-align:right;"> delta_elpd </th>
   <th style="text-align:right;"> se_delta_elpd </th>
   <th style="text-align:right;"> z </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> null </td>
   <td style="text-align:left;"> sharing </td>
   <td style="text-align:right;"> -0.423 </td>
   <td style="text-align:right;"> 1.121 </td>
   <td style="text-align:right;"> -0.377 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> null </td>
   <td style="text-align:left;"> causal </td>
   <td style="text-align:right;"> -1.968 </td>
   <td style="text-align:right;"> 2.445 </td>
   <td style="text-align:right;"> -0.805 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> sharing </td>
   <td style="text-align:left;"> causal </td>
   <td style="text-align:right;"> -1.545 </td>
   <td style="text-align:right;"> 1.531 </td>
   <td style="text-align:right;"> -1.009 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> null </td>
   <td style="text-align:left;"> sharing </td>
   <td style="text-align:right;"> 0.259 </td>
   <td style="text-align:right;"> 0.829 </td>
   <td style="text-align:right;"> 0.313 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> null </td>
   <td style="text-align:left;"> causal </td>
   <td style="text-align:right;"> -0.145 </td>
   <td style="text-align:right;"> 1.755 </td>
   <td style="text-align:right;"> -0.083 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> sharing </td>
   <td style="text-align:left;"> causal </td>
   <td style="text-align:right;"> -0.405 </td>
   <td style="text-align:right;"> 1.280 </td>
   <td style="text-align:right;"> -0.316 </td>
  </tr>
</tbody>
</table>

`````
:::
:::



::: {.cell}

```{.r .cell-code}
cause_res <- imap_dfr(cause_out, function(x, arm) {
  tidy_cause(x$fit,
             exposure = if (arm == "overlapping") cfg$exposures$glgc2021_ldl$label
                        else cfg$exposures$willer2013_ldl$label,
             outcome  = cfg$outcomes$ebmd_morris2019$label,
             arm      = arm,
             ci_size  = cfg$cause$ci_size)
})
write_csv(cause_res, pp(cfg$paths$out, "cause_results.csv"))

cause_res |>
  select(arm, b, lci, uci, eta_med, q_med, delta_elpd, se_delta_elpd, z, pval, verdict) |>
  kable(caption = "CAUSE. b is the posterior median gamma (causal model); eta is the shared-factor effect and q the proportion of variants acting through it.",
        digits = 4) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>CAUSE. b is the posterior median gamma (causal model); eta is the shared-factor effect and q the proportion of variants acting through it.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> b </th>
   <th style="text-align:right;"> lci </th>
   <th style="text-align:right;"> uci </th>
   <th style="text-align:right;"> eta_med </th>
   <th style="text-align:right;"> q_med </th>
   <th style="text-align:right;"> delta_elpd </th>
   <th style="text-align:right;"> se_delta_elpd </th>
   <th style="text-align:right;"> z </th>
   <th style="text-align:right;"> pval </th>
   <th style="text-align:left;"> verdict </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0403 </td>
   <td style="text-align:right;"> -0.0793 </td>
   <td style="text-align:right;"> -0.0026 </td>
   <td style="text-align:right;"> -0.0503 </td>
   <td style="text-align:right;"> 0.0337 </td>
   <td style="text-align:right;"> -1.5449 </td>
   <td style="text-align:right;"> 1.5314 </td>
   <td style="text-align:right;"> -1.0089 </td>
   <td style="text-align:right;"> 0.1565 </td>
   <td style="text-align:left;"> sharing not rejected </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0215 </td>
   <td style="text-align:right;"> -0.0519 </td>
   <td style="text-align:right;"> 0.0074 </td>
   <td style="text-align:right;"> -0.2020 </td>
   <td style="text-align:right;"> 0.0197 </td>
   <td style="text-align:right;"> -0.4047 </td>
   <td style="text-align:right;"> 1.2797 </td>
   <td style="text-align:right;"> -0.3162 </td>
   <td style="text-align:right;"> 0.3759 </td>
   <td style="text-align:left;"> sharing not rejected </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 2b — Bayesian reading

The ELPD z-test above is a frequentist wrapper on a model that is Bayesian
underneath. It also conflates two distinct questions. Separating them is more
informative than either p-value.

### Question 1 — given the causal model, what is $\gamma$?


::: {.cell}

```{.r .cell-code}
qtabs <- imap_dfr(cause_out, function(x, arm) {
  cause_posterior_quantiles(x$fit, model = "causal") |> mutate(arm = arm, .before = 1)
})
write_csv(qtabs, pp(cfg$paths$out, "cause_posterior_quantiles.csv"))

qtabs |>
  filter(param == "gamma") |>
  select(arm, q0.025, q0.25, q0.5, q0.75, q0.975) |>
  kable(caption = "Posterior quantiles of gamma under the causal model. The interquartile range (q0.25-q0.75) is the honest 'where is the effect' statement; no null hypothesis required.",
        digits = 4) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Posterior quantiles of gamma under the causal model. The interquartile range (q0.25-q0.75) is the honest 'where is the effect' statement; no null hypothesis required.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> q0.025 </th>
   <th style="text-align:right;"> q0.25 </th>
   <th style="text-align:right;"> q0.5 </th>
   <th style="text-align:right;"> q0.75 </th>
   <th style="text-align:right;"> q0.975 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0793 </td>
   <td style="text-align:right;"> -0.0522 </td>
   <td style="text-align:right;"> -0.0403 </td>
   <td style="text-align:right;"> -0.0279 </td>
   <td style="text-align:right;"> -0.0026 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0519 </td>
   <td style="text-align:right;"> -0.0315 </td>
   <td style="text-align:right;"> -0.0215 </td>
   <td style="text-align:right;"> -0.0120 </td>
   <td style="text-align:right;"> 0.0074 </td>
  </tr>
</tbody>
</table>

`````
:::

```{.r .cell-code}
p_neg <- imap_dfr(cause_out, function(x, arm) {
  q <- cause_posterior_quantiles(x$fit, model = "causal")
  tibble(arm = arm,
         `P(gamma < 0)`      = cause_p_below(q, "gamma", 0),
         `P(gamma < -0.02)`  = cause_p_below(q, "gamma", -0.02),
         `P(gamma < -0.05)`  = cause_p_below(q, "gamma", -0.05))
})
write_csv(p_neg, pp(cfg$paths$out, "cause_posterior_probabilities.csv"))
p_neg |>
  kable(caption = "Posterior probability that the causal effect is below a threshold, conditional on the causal model.",
        digits = 3) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Posterior probability that the causal effect is below a threshold, conditional on the causal model.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> P(gamma &lt; 0) </th>
   <th style="text-align:right;"> P(gamma &lt; -0.02) </th>
   <th style="text-align:right;"> P(gamma &lt; -0.05) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> 0.975 </td>
   <td style="text-align:right;"> 0.824 </td>
   <td style="text-align:right;"> 0.297 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> 0.929 </td>
   <td style="text-align:right;"> 0.539 </td>
   <td style="text-align:right;"> 0.032 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


### Question 2 — is the causal model the right model?

Pseudo-BMA weights spread the predictive support across all three models instead
of testing one pairwise contrast, which suits an inconclusive comparison far
better than a p-value.


::: {.cell}

```{.r .cell-code}
wts <- imap_dfr(cause_out, function(x, arm) {
  el <- as_tibble(x$fit$elpd)
  cause_model_weights(el) |> mutate(arm = arm, .before = 1)
})
write_csv(wts, pp(cfg$paths$out, "cause_model_weights.csv"))
wts |>
  select(arm, model, elpd_rel_null, weight) |>
  pivot_wider(id_cols = arm, names_from = model, values_from = weight) |>
  kable(caption = "Pseudo-BMA weights: share of predictive support for each model.",
        digits = 3) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Pseudo-BMA weights: share of predictive support for each model.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> null </th>
   <th style="text-align:right;"> sharing </th>
   <th style="text-align:right;"> causal </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> 0.103 </td>
   <td style="text-align:right;"> 0.158 </td>
   <td style="text-align:right;"> 0.739 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> 0.342 </td>
   <td style="text-align:right;"> 0.264 </td>
   <td style="text-align:right;"> 0.395 </td>
  </tr>
</tbody>
</table>

`````
:::

```{.r .cell-code}
# The weights above ignore the SE on delta_elpd, which here is comparable to
# delta_elpd itself. Propagate it so the weights carry their own uncertainty.
wts_boot <- imap_dfr(cause_out, function(x, arm) {
  cause_model_weights_boot(as_tibble(x$fit$elpd)) |> mutate(arm = arm, .before = 1)
})
write_csv(wts_boot, pp(cfg$paths$out, "cause_model_weights_boot.csv"))
wts_boot |>
  kable(caption = "Model weights with uncertainty propagated from se_delta_elpd. Wide intervals here are the real message: the data do not identify which model is right.",
        digits = 3) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Model weights with uncertainty propagated from se_delta_elpd. Wide intervals here are the real message: the data do not identify which model is right.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:left;"> model </th>
   <th style="text-align:right;"> weight_med </th>
   <th style="text-align:right;"> weight_lo </th>
   <th style="text-align:right;"> weight_hi </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> causal </td>
   <td style="text-align:right;"> 0.714 </td>
   <td style="text-align:right;"> 0.018 </td>
   <td style="text-align:right;"> 0.997 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> null </td>
   <td style="text-align:right;"> 0.085 </td>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 0.574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> sharing </td>
   <td style="text-align:right;"> 0.140 </td>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 0.829 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> causal </td>
   <td style="text-align:right;"> 0.377 </td>
   <td style="text-align:right;"> 0.018 </td>
   <td style="text-align:right;"> 0.955 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> null </td>
   <td style="text-align:right;"> 0.294 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0.703 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> sharing </td>
   <td style="text-align:right;"> 0.227 </td>
   <td style="text-align:right;"> 0.014 </td>
   <td style="text-align:right;"> 0.697 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Step 3 — posterior plots


::: {.cell}

```{.r .cell-code}
walk2(cause_out, names(cause_out), function(x, arm) {
  cat("\n\n### ", arm, "\n\n")
  print(plot(x$fit))
})
```

::: {.cell-output .cell-output-stdout}

```


###  overlapping 
```


:::

::: {.cell-output-display}
![](robust_mr_cause_files/figure-html/plots-1.png){width=864}
:::

::: {.cell-output .cell-output-stdout}

```
NULL


###  overlap-free 
```


:::

::: {.cell-output-display}
![](robust_mr_cause_files/figure-html/plots-2.png){width=864}
:::

::: {.cell-output .cell-output-stdout}

```
NULL
```


:::
:::


## Interpreting the three parameters

| Parameter | Meaning | What a large value implies |
|---|---|---|
| $\gamma$ | Causal effect of LDL-C on eBMD | The effect we are after |
| $\eta$ | Effect of the shared factor $U$ | Correlated pleiotropy is present |
| $q$ | Proportion of variants acting through $U$ | How much of the instrument set is contaminated |

The failure mode to watch for is $\gamma \approx 0$ with $\eta$ and $q$ both
substantial: that is a signal driven entirely by a shared factor, and it would mean
the conventional IVW estimate of −0.051 is confounded rather than causal.

The opposite result — sharing decisively rejected in favour of causal in **both**
arms — would be the strongest statement this project can make about
cholesterol → BMD, because it survives the one criticism the existing IVW/Egger
battery cannot answer.

A third possibility deserves saying out loud in advance: CAUSE returning "sharing
not rejected" while MR-APSS returns a clear non-zero $\beta$. That is a known
pattern — CAUSE trades power for specificity against correlated pleiotropy — and
should be reported as a discrepancy, not resolved by picking the friendlier
answer.
