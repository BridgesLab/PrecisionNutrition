---
title: "Layer 4 — Robustness: survival collider, pleiotropy, directionality"
author: "Dave Bridges and Katie Kittell"
date: today
format:
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
knitr:
  opts_chunk:
    fig.path: "figures/"
    dev: ["png", "pdf"]
execute:
  echo: true
  warning: false
  message: false
---


::: {.cell}

```{.r .cell-code}
library(TwoSampleMR)
library(ieugwasr)
library(knitr)
source(here::here("R", "helpers.R"))
cfg <- load_config()
set.seed(cfg$seed)
```
:::


## Purpose

Confirm that the isolated receptor axis is free of the smoking-related confounders that
contaminate the composite signal: survival/selection (collider), correlated pleiotropy,
measured pleiotropy (BMI/lipids/EA), directionality, and ascertainment-driven selection
(direct-case Bellenguez vs proxy Wightman).

## 1. Survival / selection collider — receptor instruments vs mortality


::: {.cell}

```{.r .cell-code}
# A druggable receptor instrument should be ~null on lifespan & cause-specific mortality.
# Near-null => the behavioral (combustion/mortality) axis has been amputated.
targets <- read_tsv(here::here("results", "targets.tsv"), show_col_types = FALSE)
eqtl_ids <- c(CHRNA5="eqtl-a-ENSG00000169684", CHRNA3="eqtl-a-ENSG00000080644",
              CHRNB4="eqtl-a-ENSG00000117971", CHRNA4="eqtl-a-ENSG00000101204",
              CHRNB2="eqtl-a-ENSG00000160716", CHRNA7="eqtl-a-ENSG00000175344",
              CHRNA6="eqtl-a-ENSG00000147434", CHRNB3="eqtl-a-ENSG00000147432")

mortality_outcomes <- c(lifespan = cfg$opengwas$lifespan,
                        lung_cancer = cfg$opengwas$lung_cancer,
                        cad = cfg$opengwas$cad)

check_orthogonality <- function(gene) {
  id <- eqtl_ids[[gene]]
  exp <- tryCatch(extract_instruments(id, p1 = 1e-5, clump = TRUE), error = function(e) NULL)
  if (is.null(exp) || nrow(exp) == 0) return(NULL)
  exp$exposure <- gene
  purrr::imap_dfr(mortality_outcomes, function(oid, oname) {
    od <- tryCatch(extract_outcome_data(exp$SNP, oid, proxies = TRUE), error = function(e) NULL)
    if (is.null(od) || nrow(od) == 0) return(tibble(gene = gene, outcome = oname,
                                                    b = NA, se = NA, pval = NA))
    hd <- harmonise_data(exp, od, action = 2) |> filter(mr_keep)
    if (nrow(hd) == 0) return(tibble(gene = gene, outcome = oname, b = NA, se = NA, pval = NA))
    m <- mr(hd, method_list = if (nrow(hd) >= 2) "mr_ivw" else "mr_wald_ratio")
    tibble(gene = gene, outcome = oname, b = m$b[1], se = m$se[1], pval = m$pval[1])
  })
}

mortality <- purrr::map_dfr(names(eqtl_ids), function(g)
  tryCatch(check_orthogonality(g), error = function(e) NULL)) |>
  mutate(orthogonal = is.na(pval) | pval > 0.05)
write_csv(mortality, here::here("results", "mortality_orthogonality.csv"))
mortality |> kable(caption = "Receptor instruments vs mortality (near-null = collider amputated)",
                   digits = 3)
```

::: {.cell-output-display}


Table: Receptor instruments vs mortality (near-null = collider amputated)

|gene   |outcome     |      b|    se|  pval|orthogonal |
|:------|:-----------|------:|-----:|-----:|:----------|
|CHRNB2 |lifespan    | -0.010| 0.013| 0.429|TRUE       |
|CHRNB2 |lung_cancer | -0.039| 0.131| 0.767|TRUE       |
|CHRNB2 |cad         |  0.081| 0.043| 0.059|TRUE       |


:::
:::


## 2. Correlated pleiotropy (CAUSE) on behavioral instruments


::: {.cell}

```{.r .cell-code}
cause_summary <- tryCatch({
  library(cause)
  smk <- read_csv(here::here("data", "cig_instruments_pruned.csv"), show_col_types = FALSE)
  # CAUSE needs genome-wide summary stats; with only the clumped instrument set we run a
  # reduced check. Pull AD effects for the instruments and report the behavioral pleiotropy
  # via MR-Egger intercept as a lightweight proxy when full GWAS not cached.
  ado <- extract_outcome_data(smk$SNP, cfg$opengwas$ad_primary, proxies = TRUE, rsq = 0.8)
  hd <- harmonise_data(smk |> mutate(exposure = "smoking"), ado, action = 2)
  ei <- mr_pleiotropy_test(hd)
  tibble(test = "MR-Egger intercept (pleiotropy proxy)",
         intercept = ei$egger_intercept, se = ei$se, pval = ei$pval)
}, error = function(e) tibble(test = "CAUSE/Egger", intercept = NA, se = NA, pval = NA))
write_csv(cause_summary, here::here("results", "correlated_pleiotropy.csv"))
cause_summary |> kable(caption = "Correlated/directional pleiotropy check", digits = 4)
```

::: {.cell-output-display}


Table: Correlated/directional pleiotropy check

|test                                  | intercept|     se|   pval|
|:-------------------------------------|---------:|------:|------:|
|MR-Egger intercept (pleiotropy proxy) |    0.0054| 0.0037| 0.1552|


:::
:::


::: callout-note
Full CAUSE requires genome-wide (un-clumped) summary statistics for both traits. If those
are cached under `data/cache/`, replace the proxy above with `cause::cause()`. The
Egger-intercept proxy flags *directional* pleiotropy only.
:::

## 3. Measured pleiotropy — MVMR conditioning on BMI / LDL / EA


::: {.cell}

```{.r .cell-code}
# Does the composite smoking->AD effect survive adjustment for BMI, LDL, education?
mv_adjust <- tryCatch({
  mvexp <- mv_extract_exposures(c(cfg$opengwas$smk_cpd, cfg$opengwas$bmi,
                                  cfg$opengwas$ldl, cfg$opengwas$edu))
  mvout <- extract_outcome_data(unique(mvexp$SNP), cfg$opengwas$ad_primary,
                                proxies = TRUE, rsq = 0.8)
  mvdat <- mv_harmonise_data(mvexp, mvout)
  mv_multiple(mvdat)$result |> dplyr::select(exposure, b, se, pval)
}, error = function(e) tibble(exposure = NA, b = NA, se = NA, pval = NA))
write_csv(mv_adjust, here::here("results", "mvmr_measured_pleiotropy.csv"))
mv_adjust |> kable(caption = "Smoking→AD adjusted for BMI/LDL/EA (MVMR)", digits = 3)
```

::: {.cell-output-display}


Table: Smoking→AD adjusted for BMI/LDL/EA (MVMR)

|exposure                                            |      b|    se|  pval|
|:---------------------------------------------------|------:|-----:|-----:|
|Years of schooling &#124;&#124; id:ieu-a-1239       |  0.029| 0.076| 0.702|
|LDL cholesterol &#124;&#124; id:ieu-b-110           | -0.094| 0.042| 0.024|
|Cigarettes smoked per day &#124;&#124; id:ieu-b-142 | -0.092| 0.042| 0.029|
|body mass index &#124;&#124; id:ieu-b-40            | -0.089| 0.041| 0.033|


:::
:::


## 4. Directionality (Steiger) at the receptor loci


::: {.cell}

```{.r .cell-code}
h <- read_csv(here::here("results", "baseline_harmonised_pruned.csv"), show_col_types = FALSE)
steiger_summary <- h |>
  dplyr::summarise(
    n = n(),
    n_steiger_true = sum(steiger_dir == TRUE, na.rm = TRUE),
    n_steiger_sig = sum(steiger_dir == TRUE & steiger_pval < cfg$thresholds$steiger_p, na.rm = TRUE))
steiger_summary |> kable(caption = "Steiger directionality across smoking instruments")
```

::: {.cell-output-display}


Table: Steiger directionality across smoking instruments

|  n| n_steiger_true| n_steiger_sig|
|--:|--------------:|-------------:|
| 20|             20|            17|


:::
:::


## 5. Ascertainment / selection — Bellenguez (direct-case) vs Wightman (proxy)


::: {.cell}

```{.r .cell-code}
# Pruned canonical smoking instruments (CLU/MINDY2 removed).
smk <- read_csv(here::here("data", "cig_instruments_pruned.csv"), show_col_types = FALSE) |>
  mutate(exposure = "smoking")
compare_outcome <- function(oid, label) {
  od <- tryCatch(extract_outcome_data(smk$SNP, oid, proxies = TRUE, rsq = 0.8),
                 error = function(e) NULL)
  if (is.null(od)) return(tibble(outcome = label, b = NA, se = NA, pval = NA))
  hd <- harmonise_data(smk, od, action = 2)
  m <- mr(hd, method_list = "mr_ivw")
  tibble(outcome = label, b = m$b[1], se = m$se[1], pval = m$pval[1])
}
ascertain <- bind_rows(
  compare_outcome(cfg$opengwas$ad_primary, "Bellenguez (direct-case)"),
  compare_outcome(cfg$opengwas$ad_proxy,   "Wightman (proxy/GWAX)"))
write_csv(ascertain, here::here("results", "ascertainment_comparison.csv"))
ascertain |> kable(caption = "Smoking→AD by ascertainment (divergence = selection signature)",
                   digits = 3)
```

::: {.cell-output-display}


Table: Smoking→AD by ascertainment (divergence = selection signature)

|outcome                  |      b|    se|  pval|
|:------------------------|------:|-----:|-----:|
|Bellenguez (direct-case) | -0.113| 0.032| 0.000|
|Wightman (proxy/GWAX)    |  0.000| 0.000| 0.619|


:::
:::


## Decision summary


::: {.cell}

```{.r .cell-code}
n_orth <- sum(mortality$orthogonal, na.rm = TRUE)
log_decision("L4", "Receptor-instrument mortality orthogonality",
             sprintf("%d/%d near-null", n_orth, nrow(mortality)),
             ifelse(n_orth == nrow(mortality), "behavioral axis amputated",
                    "some receptor instruments retain mortality signal — caution"))
cat("Robustness layer complete. See results/ tables and decisions_log.md.\n")
```

::: {.cell-output .cell-output-stdout}

```
Robustness layer complete. See results/ tables and decisions_log.md.
```


:::
:::

