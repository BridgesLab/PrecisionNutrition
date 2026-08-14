---
title: "Robust MR — Synthesis"
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

Reconcile CAUSE and MR-APSS against the conventional estimates already reported in
`ANALYSIS.md` §3, and state what the robust methods do and do not license us to
conclude about **cholesterol → BMD**.


::: {.cell}

```{.r .cell-code}
read_if <- function(f) {
  p <- pp(cfg$paths$out, f)
  if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL
}

apss   <- read_if("mrapss_results.csv")
cause  <- read_if("cause_results.csv")
bg     <- read_if("apss_background_parameters.csv")
ivw    <- read_if("ivw_baseline.csv")
raps   <- read_if("raps_results.csv")
mrbee  <- read_if("mrbee_results.csv")
rxy    <- read_if("mrbee_error_cov.csv")
```
:::


## 0. What each method assumes

The five methods split into two classes by a single assumption: whether horizontal
pleiotropic effects may be **correlated with the instrument–exposure effects**.
Methods in the first class assume they are not (the InSIDE assumption); methods in
the second explicitly model a shared confounder. Reading them together localises
where any attenuation comes from.

| Method | Uncorrelated pleiotropy | **Correlated pleiotropy** | Sample overlap | Weak instruments |
|---|---|---|---|---|
| IVW | ✗ | ✗ | ✗ | ✗ |
| MR-RAPS | ✓ random effects $\tau^2$ | ✗ | ✗ | ✓ profile score |
| MRBEE | ✓ iterative outlier removal | ✗ | ✓ error covariance $R_{xy}$ | ✓ bias-corrected |
| CAUSE | ✓ | ✓ shared factor ($\eta$, $q$) | ✓ $\rho$ | partial |
| MR-APSS | ✓ foreground | ✓ background $\Omega$ | ✓ $C$ matrix | ✓ selection correction |
| *MRAID* | *✓* | *✓ spike-slab* | *✗ — why it is excluded* | *✓* |

### Strengths and weaknesses

| Method | Strength | Weakness |
|---|---|---|
| **IVW** | Efficient, transparent, the field's default | Assumes no pleiotropy at all; SEs too narrow when any exists |
| **MR-RAPS** | Isolates the weak-instrument correction — nothing else changes. Reports $\tau^2$, so balanced pleiotropy is visible rather than assumed away | Blind to correlated pleiotropy and to overlap. Silently biased if a shared factor exists |
| **MRBEE** | Corrects overlap *and* weak instruments in one estimating equation; $R_{xy}$ is an independent overlap measurement not resting on LDSC assumptions | Pleiotropy handled by deletion, so the estimate is conditional on the outlier test being right. No shared-confounder model |
| **CAUSE** | The only method that directly tests causal vs shared-factor explanations. Conservative by design, so a positive result is meaningful | Badly underpowered at realistic effect sizes; a null is close to uninformative. Reports model comparison, not a clean effect estimate |
| **MR-APSS** | Handles all three problems at once; the relaxed 5e-5 threshold is a genuine power gain because selection bias is corrected | Inherits LDSC's assumptions. When the outcome GWAS has a large intercept (eBMD: 1.51) its SEs inflate substantially — honest, but possibly over-conservative |

## 1. Did we actually have a sample-overlap problem?

Two independent measurements of the same thing, from two different model families.
They should agree.


::: {.cell}

```{.r .cell-code}
overlap_evidence <- bind_rows(
  if (!is.null(bg)) bg |> transmute(arm, source = "MR-APSS cross-trait LDSC intercept (C12)",
                                    value = C12) else NULL,
  if (file.exists(pp(cfg$paths$out, "cause_nuisance.csv")))
    read_csv(pp(cfg$paths$out, "cause_nuisance.csv"), show_col_types = FALSE) |>
      transmute(arm, source = "CAUSE rho", value = rho) else NULL,
  if (!is.null(rxy)) rxy |> transmute(arm, source = "MRBEE error covariance (Rxy off-diagonal)",
                                      value = cov_overlap) else NULL
)

if (nrow(overlap_evidence) > 0) {
  overlap_evidence |>
    pivot_wider(names_from = arm, values_from = value) |>
    kable(caption = "Three independent estimates of shared sample structure, from three model families. All should be near zero in the overlap-free (GLGC 2013) arm and non-zero in the overlapping arm. MRBEE's estimate is the useful cross-check because, unlike the other two, it does not rest on LDSC assumptions.",
          digits = 4) |>
    kable_styling(full_width = FALSE)
}
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Three independent estimates of shared sample structure, from three model families. All should be near zero in the overlap-free (GLGC 2013) arm and non-zero in the overlapping arm. MRBEE's estimate is the useful cross-check because, unlike the other two, it does not rest on LDSC assumptions.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> source </th>
   <th style="text-align:right;"> overlapping </th>
   <th style="text-align:right;"> overlap-free </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MR-APSS cross-trait LDSC intercept (C12) </td>
   <td style="text-align:right;"> 0.0226 </td>
   <td style="text-align:right;"> -0.0074 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAUSE rho </td>
   <td style="text-align:right;"> -0.0109 </td>
   <td style="text-align:right;"> -0.0146 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MRBEE error covariance (Rxy off-diagonal) </td>
   <td style="text-align:right;"> -0.0045 </td>
   <td style="text-align:right;"> -0.0023 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## 2. All estimates side by side


::: {.cell}

```{.r .cell-code}
combined <- bind_rows(
  if (!is.null(ivw)) ivw |> transmute(method = "IVW (5e-8, this pipeline)",
                                      arm = if_else(exposure == "glgc2021_ldl",
                                                    "overlapping", "overlap-free"),
                                      b, se, lci = b - 1.96 * se, uci = b + 1.96 * se,
                                      pval, class = "no pleiotropy model") else NULL,
  if (!is.null(raps)) raps |> transmute(method = paste0(method, " @ ", threshold),
                                        arm, b, se, lci, uci, pval,
                                        class = "InSIDE (uncorrelated pleiotropy)") else NULL,
  if (!is.null(mrbee)) mrbee |> transmute(method = paste0("MRBEE @ ", threshold),
                                          arm, b, se, lci, uci, pval,
                                          class = "InSIDE (uncorrelated pleiotropy)") else NULL,
  if (!is.null(apss)) apss |> transmute(method = paste0("MR-APSS (", C_setting, ")"),
                                        arm, b, se, lci, uci, pval,
                                        class = "correlated pleiotropy modelled") else NULL,
  if (!is.null(cause)) cause |> transmute(method = "CAUSE (gamma, causal model)",
                                          arm, b, se = NA_real_, lci, uci, pval,
                                          class = "correlated pleiotropy modelled") else NULL
)
write_csv(combined, pp(cfg$paths$out, "robust_mr_combined.csv"))

combined |>
  arrange(arm, class, method) |>
  kable(caption = "LDL-C -> heel eBMD, SD per SD. The CAUSE p-value is the sharing-vs-causal ELPD test, not a Wald test on gamma - do not read it as the same quantity as the others.",
        digits = 4) |>
  kable_styling(full_width = FALSE) |>
  pack_rows(index = table(combined |> arrange(arm, class, method) |> pull(arm)))
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>LDL-C -&gt; heel eBMD, SD per SD. The CAUSE p-value is the sharing-vs-causal ELPD test, not a Wald test on gamma - do not read it as the same quantity as the others.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> method </th>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:right;"> b </th>
   <th style="text-align:right;"> se </th>
   <th style="text-align:right;"> lci </th>
   <th style="text-align:right;"> uci </th>
   <th style="text-align:right;"> pval </th>
   <th style="text-align:left;"> class </th>
  </tr>
 </thead>
<tbody>
  <tr grouplength="10"><td colspan="8" style="border-bottom: 1px solid;"><strong>overlap-free</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (l2) @ 5e-05 </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0522 </td>
   <td style="text-align:right;"> 0.0186 </td>
   <td style="text-align:right;"> -0.0887 </td>
   <td style="text-align:right;"> -0.0158 </td>
   <td style="text-align:right;"> 0.0050 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (l2) @ 5e-08 </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0587 </td>
   <td style="text-align:right;"> 0.0235 </td>
   <td style="text-align:right;"> -0.1048 </td>
   <td style="text-align:right;"> -0.0127 </td>
   <td style="text-align:right;"> 0.0123 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (tukey) @ 5e-05 </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0300 </td>
   <td style="text-align:right;"> 0.0150 </td>
   <td style="text-align:right;"> -0.0593 </td>
   <td style="text-align:right;"> -0.0006 </td>
   <td style="text-align:right;"> 0.0454 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (tukey) @ 5e-08 </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0303 </td>
   <td style="text-align:right;"> 0.0175 </td>
   <td style="text-align:right;"> -0.0645 </td>
   <td style="text-align:right;"> 0.0040 </td>
   <td style="text-align:right;"> 0.0831 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MRBEE @ 5e-05 </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0332 </td>
   <td style="text-align:right;"> 0.0164 </td>
   <td style="text-align:right;"> -0.0653 </td>
   <td style="text-align:right;"> -0.0012 </td>
   <td style="text-align:right;"> 0.0420 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MRBEE @ 5e-08 </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0366 </td>
   <td style="text-align:right;"> 0.0195 </td>
   <td style="text-align:right;"> -0.0748 </td>
   <td style="text-align:right;"> 0.0016 </td>
   <td style="text-align:right;"> 0.0607 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> CAUSE (gamma, causal model) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0215 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.0519 </td>
   <td style="text-align:right;"> 0.0074 </td>
   <td style="text-align:right;"> 0.3759 </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-APSS (corrected) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0375 </td>
   <td style="text-align:right;"> 0.0187 </td>
   <td style="text-align:right;"> -0.0741 </td>
   <td style="text-align:right;"> -0.0008 </td>
   <td style="text-align:right;"> 0.0449 </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-APSS (uncorrected) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0361 </td>
   <td style="text-align:right;"> 0.0184 </td>
   <td style="text-align:right;"> -0.0721 </td>
   <td style="text-align:right;"> -0.0001 </td>
   <td style="text-align:right;"> 0.0495 </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> IVW (5e-8, this pipeline) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0539 </td>
   <td style="text-align:right;"> 0.0068 </td>
   <td style="text-align:right;"> -0.0672 </td>
   <td style="text-align:right;"> -0.0406 </td>
   <td style="text-align:right;"> 0.0000 </td>
   <td style="text-align:left;"> no pleiotropy model </td>
  </tr>
  <tr grouplength="10"><td colspan="8" style="border-bottom: 1px solid;"><strong>overlapping</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (l2) @ 5e-05 </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0564 </td>
   <td style="text-align:right;"> 0.0209 </td>
   <td style="text-align:right;"> -0.0974 </td>
   <td style="text-align:right;"> -0.0154 </td>
   <td style="text-align:right;"> 0.0070 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (l2) @ 5e-08 </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0836 </td>
   <td style="text-align:right;"> 0.0251 </td>
   <td style="text-align:right;"> -0.1329 </td>
   <td style="text-align:right;"> -0.0343 </td>
   <td style="text-align:right;"> 0.0009 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (tukey) @ 5e-05 </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0430 </td>
   <td style="text-align:right;"> 0.0166 </td>
   <td style="text-align:right;"> -0.0754 </td>
   <td style="text-align:right;"> -0.0105 </td>
   <td style="text-align:right;"> 0.0095 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-RAPS (tukey) @ 5e-08 </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0623 </td>
   <td style="text-align:right;"> 0.0204 </td>
   <td style="text-align:right;"> -0.1023 </td>
   <td style="text-align:right;"> -0.0224 </td>
   <td style="text-align:right;"> 0.0022 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MRBEE @ 5e-05 </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0510 </td>
   <td style="text-align:right;"> 0.0195 </td>
   <td style="text-align:right;"> -0.0892 </td>
   <td style="text-align:right;"> -0.0128 </td>
   <td style="text-align:right;"> 0.0089 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MRBEE @ 5e-08 </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0715 </td>
   <td style="text-align:right;"> 0.0216 </td>
   <td style="text-align:right;"> -0.1137 </td>
   <td style="text-align:right;"> -0.0292 </td>
   <td style="text-align:right;"> 0.0009 </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> CAUSE (gamma, causal model) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0403 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.0793 </td>
   <td style="text-align:right;"> -0.0026 </td>
   <td style="text-align:right;"> 0.1565 </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-APSS (corrected) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0561 </td>
   <td style="text-align:right;"> 0.0265 </td>
   <td style="text-align:right;"> -0.1080 </td>
   <td style="text-align:right;"> -0.0041 </td>
   <td style="text-align:right;"> 0.0344 </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MR-APSS (uncorrected) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0517 </td>
   <td style="text-align:right;"> 0.0262 </td>
   <td style="text-align:right;"> -0.1030 </td>
   <td style="text-align:right;"> -0.0004 </td>
   <td style="text-align:right;"> 0.0481 </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> IVW (5e-8, this pipeline) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0672 </td>
   <td style="text-align:right;"> 0.0056 </td>
   <td style="text-align:right;"> -0.0782 </td>
   <td style="text-align:right;"> -0.0563 </td>
   <td style="text-align:right;"> 0.0000 </td>
   <td style="text-align:left;"> no pleiotropy model </td>
  </tr>
</tbody>
</table>

`````
:::
:::



::: {.cell}

```{.r .cell-code}
if (nrow(combined) > 0) {
  robust_mr_forest(combined,
                   title = "LDL-C -> heel eBMD across estimators and exposure arms") +
    ggplot2::facet_grid(class ~ arm, scales = "free_y", space = "free_y") +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
}
```

::: {.cell-output-display}
![](robust_mr_summary_files/figure-html/forest-1.png){width=864}
:::
:::


## 2b. Does the attenuation come from pleiotropy or from instruments?

This is what adding the InSIDE-class methods buys. IVW → MR-APSS/CAUSE confounds
three corrections at once. MR-RAPS and MRBEE hold the correlated-pleiotropy
assumption fixed and vary only the weak-instrument (RAPS) and
weak-instrument-plus-overlap (MRBEE) handling, so the comparison separates them.


::: {.cell}

```{.r .cell-code}
if (nrow(combined) > 0) {
  combined |>
    group_by(arm, class) |>
    summarise(n_fits = n(), b_median = median(b, na.rm = TRUE),
              b_min = min(b, na.rm = TRUE), b_max = max(b, na.rm = TRUE),
              .groups = "drop") |>
    kable(caption = "Estimates grouped by assumption class. If the InSIDE class sits close to IVW while the correlated-pleiotropy class is attenuated, the shrinkage is specifically attributable to modelling a shared factor — not to weak instruments or overlap.",
          digits = 4) |>
    kable_styling(full_width = FALSE)
}
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Estimates grouped by assumption class. If the InSIDE class sits close to IVW while the correlated-pleiotropy class is attenuated, the shrinkage is specifically attributable to modelling a shared factor — not to weak instruments or overlap.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> arm </th>
   <th style="text-align:left;"> class </th>
   <th style="text-align:right;"> n_fits </th>
   <th style="text-align:right;"> b_median </th>
   <th style="text-align:right;"> b_min </th>
   <th style="text-align:right;"> b_max </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> -0.0349 </td>
   <td style="text-align:right;"> -0.0587 </td>
   <td style="text-align:right;"> -0.0300 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> -0.0361 </td>
   <td style="text-align:right;"> -0.0375 </td>
   <td style="text-align:right;"> -0.0215 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:left;"> no pleiotropy model </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0539 </td>
   <td style="text-align:right;"> -0.0539 </td>
   <td style="text-align:right;"> -0.0539 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> InSIDE (uncorrelated pleiotropy) </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> -0.0594 </td>
   <td style="text-align:right;"> -0.0836 </td>
   <td style="text-align:right;"> -0.0430 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> correlated pleiotropy modelled </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> -0.0517 </td>
   <td style="text-align:right;"> -0.0561 </td>
   <td style="text-align:right;"> -0.0403 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:left;"> no pleiotropy model </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0672 </td>
   <td style="text-align:right;"> -0.0672 </td>
   <td style="text-align:right;"> -0.0672 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


Decision rule, fixed in advance:

- **InSIDE class ≈ IVW, correlated class attenuated** → the attenuation is the
  shared-factor modelling. Correlated pleiotropy is the live concern, and CAUSE's
  inability to exclude it is the binding limitation.
- **InSIDE class also attenuated** → the shrinkage is about instrument quality
  (weak instruments, winner's curse), and says nothing about pleiotropy. The
  correlated-pleiotropy methods are then not adding what they appear to add.
- **MRBEE ≠ MR-RAPS** → the difference is overlap plus outlier deletion, since
  that is all that separates them. Cross-check against the $C_{12}$/$\rho$/$R_{xy}$
  table above.

## 3. Reading the result grid

The four informative patterns, decided in advance so the interpretation is not
chosen after seeing the numbers:

| MR-APSS | CAUSE | Reading |
|---|---|---|
| $\beta<0$, significant | causal preferred | The strongest available support for a genuine cholesterol → BMD effect. Survives overlap, correlated pleiotropy, and winner's curse. |
| $\beta<0$, significant | sharing not rejected | Effect is real *or* driven by a shared factor; CAUSE cannot separate them at this power. Report both; do not claim causality. |
| $\beta \approx 0$ | sharing not rejected | The conventional IVW estimate (−0.051) was substantially overlap- and/or pleiotropy-driven. This would be a material revision to `ANALYSIS.md` §3. |
| $\beta \approx 0$ | causal preferred | Unusual. Check the $C$ matrix — a large LDSC intercept inflates MR-APSS standard errors and can null out a real effect. |

Note that none of these bear on the **HMGCR cis** result (β ≈ −0.115). That is a
drug-target design with a different identifying assumption, and it is not
adjudicated by genome-wide polygenic methods. The relationship between the two is
the interesting part: if genome-wide LDL → BMD collapses under CAUSE/MR-APSS while
the HMGCR cis effect stands, that *strengthens* the paper's central claim that the
mechanism is mevalonate-specific rather than LDL-mediated.

## 4. Limitations to carry into the write-up

- **Two of the five methods assume InSIDE.** MR-RAPS and MRBEE assume horizontal
  pleiotropic effects are uncorrelated with the instrument–exposure effects. If a
  shared factor drives both LDL-C and eBMD, both are biased and neither will warn
  you — that failure mode is invisible to them by construction. They are included
  precisely because their assumptions differ from CAUSE's and MR-APSS's, not
  because they are more trustworthy.
- **MRBEE removes pleiotropic instruments rather than modelling them**, so its
  estimate is conditional on the IMRP outlier test making the right calls. Compare
  against MR-RAPS with the `tukey` loss, which downweights instead of deleting.
- **MRAID was not run.** It is a strict two-sample method and the primary design
  shares UK Biobank across exposure and outcome. It is deferred to a UK Biobank →
  MGI/BioVU design, mirroring the architecture already used for serum calcium.
- **Both methods assume the LDSC/polygenic model.** MR-APSS inherits LDSC's
  assumptions directly; CAUSE assumes the shared factor acts on both traits
  through the same variants. Neither handles a *reverse* effect of BMD on lipids —
  run the reverse direction separately if that is in question.
- **European ancestry only.** `eur_w_ld_chr` and the 1000G EUR panel restrict this
  to the European subsets; the multi-ancestry power of GLGC 2021 is not used.
- **CAUSE is conservative by design.** A null CAUSE result is weak evidence of no
  effect, not strong evidence.
- **The eBMD outcome is heel quantitative ultrasound**, not DXA. The femoral-neck
  DXA result in `ANALYSIS.md` was null but underpowered (GEFOS n ≈ 33k), and that
  cohort is too small for LDSC-based methods.

## 5. Suggested text for ANALYSIS.md

```
### 3b. Robustness of cholesterol -> BMD to overlap and correlated pleiotropy
_Scripts: robust_mr_prep.qmd, robust_mr_apss.qmd, robust_mr_cause.qmd._

The heel BMD result rests on UK Biobank on both sides of the analysis. Two
genome-wide methods were applied that model this explicitly: MR-APSS (sample
structure via bivariate LDSC intercepts, plus a winner's-curse correction that
licenses a relaxed 5e-5 instrument threshold) and CAUSE (shared-factor model for
correlated horizontal pleiotropy, robust to overlap via rho). Each was run on an
overlapping exposure (GLGC 2021) and an overlap-free exposure (GLGC 2013,
pre-UK-Biobank).

[FILL: C12 = _, rho = _ in the overlapping arm vs _ and _ in the overlap-free arm]
[FILL: MR-APSS beta = _ (95% CI _, _), vs _ with the correction disabled]
[FILL: CAUSE delta_elpd = _, p = _ for sharing vs causal]

MRAID was considered and excluded: it assumes non-overlapping samples, which this
design violates by construction. It is planned for a UK Biobank -> MGI/BioVU
replication.
```
