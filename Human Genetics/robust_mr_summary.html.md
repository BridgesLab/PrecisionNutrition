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
```
:::


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
      transmute(arm, source = "CAUSE rho", value = rho) else NULL
)

if (nrow(overlap_evidence) > 0) {
  overlap_evidence |>
    pivot_wider(names_from = arm, values_from = value) |>
    kable(caption = "Two estimates of shared sample structure. Both should be near zero in the overlap-free (GLGC 2013) arm and non-zero in the overlapping (GLGC 2021 + UKB eBMD) arm. If C12 and rho disagree in sign or magnitude, stop and check the harmonisation.",
          digits = 4) |>
    kable_styling(full_width = FALSE)
}
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Two estimates of shared sample structure. Both should be near zero in the overlap-free (GLGC 2013) arm and non-zero in the overlapping (GLGC 2021 + UKB eBMD) arm. If C12 and rho disagree in sign or magnitude, stop and check the harmonisation.</caption>
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
                                      pval) else NULL,
  if (!is.null(apss)) apss |> transmute(method = paste0("MR-APSS (", C_setting, ")"),
                                        arm, b, se, lci, uci, pval) else NULL,
  if (!is.null(cause)) cause |> transmute(method = "CAUSE (gamma, causal model)",
                                          arm, b, se = NA_real_, lci, uci, pval) else NULL
)
write_csv(combined, pp(cfg$paths$out, "robust_mr_combined.csv"))

combined |>
  arrange(arm, method) |>
  kable(caption = "LDL-C -> heel eBMD, SD per SD. The CAUSE p-value is the sharing-vs-causal ELPD test, not a Wald test on gamma - do not read it as the same quantity as the others.",
        digits = 4) |>
  kable_styling(full_width = FALSE)
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
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CAUSE (gamma, causal model) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0215 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.0519 </td>
   <td style="text-align:right;"> 0.0074 </td>
   <td style="text-align:right;"> 0.3759 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IVW (5e-8, this pipeline) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0539 </td>
   <td style="text-align:right;"> 0.0068 </td>
   <td style="text-align:right;"> -0.0672 </td>
   <td style="text-align:right;"> -0.0406 </td>
   <td style="text-align:right;"> 0.0000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MR-APSS (corrected) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0375 </td>
   <td style="text-align:right;"> 0.0187 </td>
   <td style="text-align:right;"> -0.0741 </td>
   <td style="text-align:right;"> -0.0008 </td>
   <td style="text-align:right;"> 0.0449 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MR-APSS (uncorrected) </td>
   <td style="text-align:left;"> overlap-free </td>
   <td style="text-align:right;"> -0.0361 </td>
   <td style="text-align:right;"> 0.0184 </td>
   <td style="text-align:right;"> -0.0721 </td>
   <td style="text-align:right;"> -0.0001 </td>
   <td style="text-align:right;"> 0.0495 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAUSE (gamma, causal model) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0403 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.0793 </td>
   <td style="text-align:right;"> -0.0026 </td>
   <td style="text-align:right;"> 0.1565 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IVW (5e-8, this pipeline) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0672 </td>
   <td style="text-align:right;"> 0.0056 </td>
   <td style="text-align:right;"> -0.0782 </td>
   <td style="text-align:right;"> -0.0563 </td>
   <td style="text-align:right;"> 0.0000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MR-APSS (corrected) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0561 </td>
   <td style="text-align:right;"> 0.0265 </td>
   <td style="text-align:right;"> -0.1080 </td>
   <td style="text-align:right;"> -0.0041 </td>
   <td style="text-align:right;"> 0.0344 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MR-APSS (uncorrected) </td>
   <td style="text-align:left;"> overlapping </td>
   <td style="text-align:right;"> -0.0517 </td>
   <td style="text-align:right;"> 0.0262 </td>
   <td style="text-align:right;"> -0.1030 </td>
   <td style="text-align:right;"> -0.0004 </td>
   <td style="text-align:right;"> 0.0481 </td>
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
                   title = "LDL-C -> heel eBMD across estimators and exposure arms")
}
```

::: {.cell-output-display}
![](robust_mr_summary_files/figure-html/forest-1.png){width=768}
:::
:::


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
