---
title: "Robust MR — MR-RAPS (uncorrelated pleiotropy, weak instruments)"
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

MR-RAPS (Zhao et al., *Ann. Statist.* 2020) sits deliberately between IVW and the
confounder-modelling methods. It assumes horizontal pleiotropy exists and may be
large, models it as random effects $\alpha_j \sim N(0, \tau^2)$, and — critically —
assumes those effects are **independent of the instrument–exposure effects**
(the InSIDE assumption). It does *not* posit a shared confounder.

Two properties earn it a place here:

1. **Weak-instrument correction by profile score.** This is RAPS's core
   contribution and the reason it tolerates a relaxed p-value threshold. It lets
   us ask whether the attenuation seen in MR-APSS is really about pleiotropy or
   just about instrument quality.
2. **$\hat\tau^2$ is reported.** Unlike IVW, RAPS estimates how much pleiotropy
   there is. A large $\tau^2$ with a stable $\beta$ is the signature of *balanced*
   pleiotropy — exactly the situation IVW handles badly and RAPS handles well.

::: callout-note
## Where this sits among the five methods

| | uncorrelated pleiotropy | **correlated** pleiotropy | overlap | weak IV |
|---|---|---|---|---|
| IVW | ✗ | ✗ | ✗ | ✗ |
| **MR-RAPS** | ✓ $\tau^2$ | ✗ | ✗ | ✓ |
| **MRBEE** | ✓ outlier removal | ✗ | ✓ | ✓ |
| CAUSE | ✓ | ✓ | ✓ | partial |
| MR-APSS | ✓ | ✓ | ✓ | ✓ |

RAPS is the only method here that corrects weak instruments *without* also
touching overlap or correlated pleiotropy. That isolation is the point: it tells
us what the weak-instrument correction alone is worth.
:::

## Data

Instruments come from the same harmonised substrate as MR-APSS and CAUSE
(`est_paras()$dat`), so differences between methods are attributable to the
method rather than to QC or harmonisation.


::: {.cell cache.extra='["d069e35e7b64ce8debe399bf833eadde","0c3b5d3f6f5ae48656f746676715cb95"]'}

```{.r .cell-code}
paras_path <- pp(cfg$paths$cache, "apss_paras.rds")
if (!file.exists(paras_path)) stop("Run robust_mr_apss.qmd first — missing ", paras_path)
paras <- readRDS(paras_path)

exp_of <- function(arm) if (arm == "overlapping") "glgc2021_ldl" else "willer2013_ldl"

map_dfr(paras, ~ tibble(n_snp = nrow(.x$dat)), .id = "arm") |>
  kable(caption = "Harmonised variants available per arm",
        format.args = list(big.mark = ",")) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Harmonised variants available per arm</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> n_snp </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> 1,017,040 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> 901,711 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Instruments at two thresholds

RAPS is run at both the conventional 5e-8 and the relaxed 5e-5 used by MR-APSS.
If the estimate is stable across thresholds, weak-instrument bias is not driving
anything.


::: {.cell cache.extra='["d069e35e7b64ce8debe399bf833eadde","0c3b5d3f6f5ae48656f746676715cb95"]'}

```{.r .cell-code}
thresholds <- c(`5e-8` = 5e-8, `5e-5` = cfg$mrapss$iv_threshold)

ivs <- purrr::map(names(paras), function(arm) {
  purrr::map(names(thresholds), function(th) {
    robust_iv_set(paras[[arm]]$dat, thresholds[[th]], cfg,
                  r2 = cfg$mrapss$clump_r2, kb = cfg$mrapss$clump_kb)
  }) |> setNames(names(thresholds))
}) |> setNames(names(paras))

map_dfr(names(ivs), function(arm) {
  map_dfr(names(thresholds), function(th) {
    d <- ivs[[arm]][[th]]
    tibble(arm = arm, threshold = th, n_iv = nrow(d),
           mean_F = mean((d$b.exp / d$se.exp)^2))
  })
}) |>
  kable(caption = "Instruments by arm and threshold. Mean F below ~10 at the relaxed threshold is precisely the regime RAPS is built for.",
        digits = 1) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Instruments by arm and threshold. Mean F below ~10 at the relaxed threshold is precisely the regime RAPS is built for.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:left;"> threshold </th>
   <th style="text-align:right;"> n_iv </th>
   <th style="text-align:right;"> mean_F </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> 5e-8 </td>
   <td style="text-align:right;"> 443 </td>
   <td style="text-align:right;"> 128.6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> 5e-5 </td>
   <td style="text-align:right;"> 948 </td>
   <td style="text-align:right;"> 71.5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> 5e-8 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 60.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> 5e-5 </td>
   <td style="text-align:right;"> 164 </td>
   <td style="text-align:right;"> 39.4 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Fit

Two loss functions: plain `l2` (efficient if there are no outliers) and `tukey`
(robust, downweights idiosyncratic pleiotropic outliers). `over.dispersion = TRUE`
estimates $\tau^2$ rather than assuming no pleiotropy.


::: {.cell cache.extra='["d069e35e7b64ce8debe399bf833eadde","0c3b5d3f6f5ae48656f746676715cb95"]'}

```{.r .cell-code}
library(mr.raps)

fit_one <- function(arm, th, loss) {
  d <- ivs[[arm]][[th]]
  if (nrow(d) < 4) return(NULL)
  res <- try(fit_raps(b_exp = d$b.exp, b_out = d$b.out,
                      se_exp = d$se.exp, se_out = d$se.out,
                      loss = loss, over.dispersion = TRUE),
             silent = TRUE)
  if (inherits(res, "try-error")) {
    message("RAPS failed for ", arm, " / ", th, " / ", loss, ":\n",
            conditionMessage(attr(res, "condition")))
    return(NULL)
  }
  tidy_raps(res, arm = arm, threshold = th, loss = loss, n_iv = nrow(d))
}

raps_res <- tidyr::expand_grid(arm = names(ivs), th = names(thresholds),
                               loss = c("l2", "tukey")) |>
  purrr::pmap_dfr(function(arm, th, loss) fit_one(arm, th, loss))
```

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-1.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.0836, standard error: 0.0251, p-value: 0.000883.
Estimated overdispersion variance: 2.45e-05, standard error: 1.8e-06, p-value: 5.13e-42.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                       Df Sum Sq Mean Sq
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  25   29.7 1.18810
Residuals                                             418  413.3 0.98875
                                                      F value Pr(>F)
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  1.2016 0.2319
Residuals                                                           
```


:::

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-2.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.0623, standard error: 0.0204, p-value: 0.00222.
Estimated overdispersion variance: 1.44e-05, standard error: 1.21e-06, p-value: 7.73e-33.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                       Df Sum Sq Mean Sq
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  25  45.17  1.8070
Residuals                                             418 666.62  1.5948
                                                      F value Pr(>F)
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  1.1331 0.3008
Residuals                                                           
```


:::

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-3.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.0564, standard error: 0.0209, p-value: 0.00702.
Estimated overdispersion variance: 1.95e-05, standard error: 1e-06, p-value: 3.89e-84.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                       Df Sum Sq Mean Sq
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  50     67 1.33999
Residuals                                             898    881 0.98107
                                                      F value  Pr(>F)  
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  1.3658 0.04954 *
Residuals                                                              
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-4.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.043, standard error: 0.0166, p-value: 0.00952.
Estimated overdispersion variance: 1.07e-05, standard error: 6.41e-07, p-value: 3.09e-62.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                       Df  Sum Sq Mean Sq
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  50  101.51  2.0302
Residuals                                             898 1490.35  1.6596
                                                      F value Pr(>F)
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  1.2233 0.1421
Residuals                                                           
```


:::

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-5.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.0587, standard error: 0.0235, p-value: 0.0123.
Estimated overdispersion variance: 2.45e-05, standard error: 4.38e-06, p-value: 2.4e-08.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                      Df Sum Sq Mean Sq F value
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  7  6.094 0.87061  0.8592
Residuals                                             68 68.906 1.01332        
                                                      Pr(>F)
bs(weights, knots = quantile(weights, 1:df/(df + 1))) 0.5431
Residuals                                                   
diagonal element is zero 
[1] 2
diagonal element is zero 
[1] 2
diagonal element is zero 
[1] 2
```


:::

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-6.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.0303, standard error: 0.0175, p-value: 0.0831.
Estimated overdispersion variance: 1.17e-05, standard error: 2.47e-06, p-value: 1.96e-06.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                      Df  Sum Sq Mean Sq
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  7  15.195  2.1707
Residuals                                             68 130.233  1.9152
                                                      F value Pr(>F)
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  1.1334 0.3529
Residuals                                                           
```


:::

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-7.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.0522, standard error: 0.0186, p-value: 0.00496.
Estimated overdispersion variance: 2.12e-05, standard error: 2.6e-06, p-value: 3.82e-16.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                       Df  Sum Sq Mean Sq
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  11   9.357 0.85062
Residuals                                             153 154.643 1.01074
                                                      F value Pr(>F)
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  0.8416 0.5989
Residuals                                                           
```


:::

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/fit-8.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```
Estimated causal effect: -0.03, standard error: 0.015, p-value: 0.0454.
Estimated overdispersion variance: 1.22e-05, standard error: 1.72e-06, p-value: 1.45e-12.
ANOVA test: are the weights and residuals independent? 
Analysis of Variance Table

Response: std.resids
                                                       Df  Sum Sq Mean Sq
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  11  17.717  1.6107
Residuals                                             153 250.407  1.6366
                                                      F value Pr(>F)
bs(weights, knots = quantile(weights, 1:df/(df + 1)))  0.9841 0.4633
Residuals                                                           
```


:::

```{.r .cell-code}
if (nrow(raps_res) == 0) {
  stop("All MR-RAPS fits failed — the messages above list what each interface reported.")
}
message("MR-RAPS interface used: ", paste(unique(raps_res$interface), collapse = ", "))
write_csv(raps_res, pp(cfg$paths$out, "raps_results.csv"))

raps_res |>
  select(arm, threshold, method, n_iv, b, se, lci, uci, pval, tau2, interface) |>
  arrange(arm, threshold, method) |>
  kable(caption = "MR-RAPS. tau2 is the estimated variance of horizontal pleiotropic effects; a large value with a stable beta indicates balanced pleiotropy.",
        digits = 4) |>
  kable_styling(full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>MR-RAPS. tau2 is the estimated variance of horizontal pleiotropic effects; a large value with a stable beta indicates balanced pleiotropy.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:left;"> threshold </th>
   <th style="text-align:left;"> method </th>
   <th style="text-align:right;"> n_iv </th>
   <th style="text-align:right;"> b </th>
   <th style="text-align:right;"> se </th>
   <th style="text-align:right;"> lci </th>
   <th style="text-align:right;"> uci </th>
   <th style="text-align:right;"> pval </th>
   <th style="text-align:right;"> tau2 </th>
   <th style="text-align:left;"> interface </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> 5e-5 </td>
   <td style="text-align:left;"> MR-RAPS (l2) </td>
   <td style="text-align:right;"> 164 </td>
   <td style="text-align:right;"> -0.0522 </td>
   <td style="text-align:right;"> 0.0186 </td>
   <td style="text-align:right;"> -0.0887 </td>
   <td style="text-align:right;"> -0.0158 </td>
   <td style="text-align:right;"> 0.0050 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> 5e-5 </td>
   <td style="text-align:left;"> MR-RAPS (tukey) </td>
   <td style="text-align:right;"> 164 </td>
   <td style="text-align:right;"> -0.0300 </td>
   <td style="text-align:right;"> 0.0150 </td>
   <td style="text-align:right;"> -0.0593 </td>
   <td style="text-align:right;"> -0.0006 </td>
   <td style="text-align:right;"> 0.0454 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> 5e-8 </td>
   <td style="text-align:left;"> MR-RAPS (l2) </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> -0.0587 </td>
   <td style="text-align:right;"> 0.0235 </td>
   <td style="text-align:right;"> -0.1048 </td>
   <td style="text-align:right;"> -0.0127 </td>
   <td style="text-align:right;"> 0.0123 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> 5e-8 </td>
   <td style="text-align:left;"> MR-RAPS (tukey) </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> -0.0303 </td>
   <td style="text-align:right;"> 0.0175 </td>
   <td style="text-align:right;"> -0.0645 </td>
   <td style="text-align:right;"> 0.0040 </td>
   <td style="text-align:right;"> 0.0831 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> 5e-5 </td>
   <td style="text-align:left;"> MR-RAPS (l2) </td>
   <td style="text-align:right;"> 948 </td>
   <td style="text-align:right;"> -0.0564 </td>
   <td style="text-align:right;"> 0.0209 </td>
   <td style="text-align:right;"> -0.0974 </td>
   <td style="text-align:right;"> -0.0154 </td>
   <td style="text-align:right;"> 0.0070 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> 5e-5 </td>
   <td style="text-align:left;"> MR-RAPS (tukey) </td>
   <td style="text-align:right;"> 948 </td>
   <td style="text-align:right;"> -0.0430 </td>
   <td style="text-align:right;"> 0.0166 </td>
   <td style="text-align:right;"> -0.0754 </td>
   <td style="text-align:right;"> -0.0105 </td>
   <td style="text-align:right;"> 0.0095 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> 5e-8 </td>
   <td style="text-align:left;"> MR-RAPS (l2) </td>
   <td style="text-align:right;"> 443 </td>
   <td style="text-align:right;"> -0.0836 </td>
   <td style="text-align:right;"> 0.0251 </td>
   <td style="text-align:right;"> -0.1329 </td>
   <td style="text-align:right;"> -0.0343 </td>
   <td style="text-align:right;"> 0.0009 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> 5e-8 </td>
   <td style="text-align:left;"> MR-RAPS (tukey) </td>
   <td style="text-align:right;"> 443 </td>
   <td style="text-align:right;"> -0.0623 </td>
   <td style="text-align:right;"> 0.0204 </td>
   <td style="text-align:right;"> -0.1023 </td>
   <td style="text-align:right;"> -0.0224 </td>
   <td style="text-align:right;"> 0.0022 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> mr.raps(data.frame) </td>
  </tr>
</tbody>
</table>

`````
:::
:::



::: {.cell}

```{.r .cell-code}
if (nrow(raps_res) > 0) {
  raps_res |>
    mutate(label = paste0(method, " @ ", threshold)) |>
    robust_mr_forest(title = "MR-RAPS: LDL-C -> heel eBMD, by threshold and loss") +
    facet_wrap(~arm)
}
```

::: {.cell-output-display}
![](robust_mr_raps_files/figure-html/plot-1.png){width=768}
:::
:::


## Reading it

- **Stable across thresholds** → weak-instrument bias is not the story, and the
  relaxed threshold in MR-APSS is safe.
- **Larger $|\beta|$ at 5e-8 than 5e-5** → the extra weak instruments carry
  proportionally more pleiotropy than signal.
- **$\tau^2$ substantially above zero** → real horizontal pleiotropy, but *balanced*.
  RAPS accommodates it; IVW's standard errors would be too narrow.
- **`tukey` and `l2` disagree** → a few outlying variants are driving the estimate;
  prefer the robust fit and check which SNPs.

What RAPS *cannot* tell you is whether pleiotropy is correlated with instrument
strength. If a shared factor drives both LDL-C and eBMD, RAPS is biased and will
not warn you. That is what CAUSE and MR-APSS are for, and why the two classes are
reported side by side in `robust_mr_summary.qmd`.
