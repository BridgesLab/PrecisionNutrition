---
title: "Step 2 — Naive exposure → AD, then SlopeHunter selection correction"
subtitle: "Adjusting each exposure → AD estimate for the survival collider flagged in Step 1"
author: "Dave Bridges and Katie Kittell"
date: today
format:
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
  pdf: default
knitr:
  opts_chunk:
    fig.path: "figures/"
    dev: ["png", "pdf"]
    fig.keep: "all"
    autodep: true
execute:
  echo: true
  warning: false
  message: false
  cache: true
---

## Purpose

Step 1 (`exposure_survival_screen.qmd`) established two things: (i) AD shortens
survival, so the survival collider `S` is real, and (ii) for the exposures
flagged there, the `E → S` arm is present — i.e. selecting on survival induces a
spurious `G_exposure`–AD association, in a direction given by
`sign(b_E · b_AD)`. This notebook does the correction: for each exposure we

1. estimate the **naive** `E → AD` effect (standard two-sample MR, *uncorrected*
   — this is the estimate that carries the survival-collider bias), then
2. subtract the selection-induced component using a **SlopeHunter** selection
   slope `b_SH`, exactly as in the smoking→AD analysis
   (`smoking/04e_overlap_indexevent.qmd`).

Because IVW is linear, the correction is exact:
$$ \text{corrected}(E\to\text{AD}) \;=\; (E\to\text{AD}) \;-\; b_{SH}\,(E\to\text{lifespan}). $$

`b_SH` is the SlopeHunter slope relating SNP effects on AD (prognosis) to SNP
effects on lifespan (the selection axis). It is a property of the **(AD GWAS ×
lifespan GWAS)** pair and is *exposure-independent* — estimate it once per AD
GWAS (`R/fit_selection_slope.R`, genome-wide, APOE excluded), then apply it to
every exposure. `E → lifespan` is reused from Step 1.

## Setup


::: {.cell}

```{.r .cell-code}
# cache=FALSE is required: a cached setup chunk is restored WITHOUT re-running library()/source(),
# so on re-render dplyr etc. would not be attached ("could not find function transmute").
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
# Optional, needed only for those naive models: MRPRESSO (MR-PRESSO), mr.raps (MR-RAPS).
# install.packages(c("MRPRESSO", "mr.raps"))  # both used via TwoSampleMR wrappers; absent -> IVW-MRE fallback

# Correction engine (correct_estimate + the genome-wide fit driver).
source("R/fit_selection_slope.R")

stopifnot(nchar(ieugwasr::get_opengwas_jwt()) > 0)
```
:::


## Configuration 


::: {.cell}

```{.r .cell-code}
# ---- Exposures + the SOP-validated naive model used for each E -> AD association. ----
# One row per exposure: OpenGWAS `id` AND the `naive_model` you validated for that association,
# so the collider correction runs on the exact model you used — not a re-derived default.
# Model choice can be a judgement call (e.g. LDL-C: visual APOE outliers -> MR-PRESSO), which is
# why it is recorded explicitly here rather than auto-selected. `note` documents the rationale.
# Allowed `naive_model`: "IVW-FE", "IVW-MRE", "MR-Egger", "Weighted median", "Weighted mode",
# "MR-RAPS", "MR-PRESSO"  (SOP: https://bridgeslab.sph.umich.edu/protocols/index.php/Mendelian_Randomization)
exposures <- tibble::tribble(
  ~label,                        ~id,          ~naive_model,  ~note,
  "Smoking (GSCAN, ieu-b-142)",  "ieu-b-142",  "MR-PRESSO",     "visual APOE outliers -> PRESSO",
  "LDL-C",                       "ieu-b-110",  "MR-PRESSO",   "visual APOE outliers -> PRESSO"
)

# ---- AD outcomes (the downstream trait, measured only in survivors). ----
ad_gwas <- c(
  "AD — Bellenguez 2022 (proxy-incl.)" = "ebi-a-GCST90027158",
  "AD — Kunkle 2019 IGAP (clinical)"   = "ieu-b-2"
)

# ---- Selection axis = the lifespan GWAS b_SH is fit on. MUST match Step 1's PRIMARY_SURVIVAL
#      and the E -> lifespan estimates we reuse below. ----
SELECTION_LABEL <- "Parental lifespan (combined, Martingale resid)"  # ebi-a-GCST006697

# ---- Selection slope b_SH per AD GWAS. ----
# b_SH is genome-wide, LD-clumped, APOE-excluded (see R/fit_selection_slope.R). It is heavy to
# fit (needs full sumstats), so we consume a cached fits CSV per AD GWAS. Bellenguez x this same
# lifespan is already fit in the smoking pipeline; Kunkle needs its own fit (chunk below).
bsh_sources <- tibble::tribble(
  ~ad_id,               ~fits_csv,                                             ~fit_set,
  "ebi-a-GCST90027158", "smoking/results/indexevent_slopehunter_fits.csv",     "excl_APOE (PRIMARY)",
  "ieu-b-2",            "R/fits/kunkle_lifespan_slopehunter_fits.csv",          "excl_APOE (PRIMARY)"
)

# Instrument / clumping parameters (match Step 1).
P_THRESH <- 5e-8
R2       <- 0.001
KB       <- 10000

# ROPE for judging whether a *corrected* effect is practically non-null (per-SD outcome).
ROPE <- 0.01
```
:::


## Naive exposure → AD (uncorrected)

The bias-carrying estimate. Same instrument selection as Step 1, pointed at AD.


::: {.cell}

```{.r .cell-code}
# Retry a (possibly rate-limited) OpenGWAS call with exponential backoff. `fn` is a thunk.
og_retry <- function(fn, tries = 4, base = 3) {
  for (i in seq_len(tries)) {
    r <- tryCatch(fn(), error = function(e) NULL)
    if (!is.null(r) && (!is.data.frame(r) || nrow(r) > 0)) return(r)
    if (i < tries) Sys.sleep(base * i)
  }
  NULL
}

# Instruments are pulled ONCE per exposure and reused across AD outcomes (halves API calls and
# keeps the instrument set identical across the two AD GWAS).
.inst_cache <- new.env(parent = emptyenv())
get_instruments <- function(exp_id, p, r2, kb) {
  key <- paste(exp_id, p, r2, kb)
  if (is.null(.inst_cache[[key]]))
    .inst_cache[[key]] <- og_retry(function() extract_instruments(exp_id, p1 = p, r2 = r2, kb = kb))
  .inst_cache[[key]]
}

# Dispatch to the SOP-validated naive model on a harmonised dataset. Returns one row with the
# estimate + diagnostics (Cochran's Q p, Egger-intercept p, and — for MR-PRESSO — #outliers and
# the global-test p). Unknown/failed methods and missing optional packages (MRPRESSO, mr.raps)
# degrade to IVW-MRE, and the `method` column shows what actually ran so a fallback is visible.
run_naive_model <- function(dat, model = "IVW-MRE") {
  tsm <- c("IVW-FE" = "mr_ivw_fe", "IVW-MRE" = "mr_ivw_mre", "MR-Egger" = "mr_egger_regression",
           "Weighted median" = "mr_weighted_median", "Weighted mode" = "mr_weighted_mode",
           "MR-RAPS" = "mr_raps")
  Q_p   <- tryCatch(mr_heterogeneity(dat, method_list = "mr_ivw")$Q_pval[1], error = \(e) NA_real_)
  egg_p <- tryCatch(mr_pleiotropy_test(dat)$pval,                            error = \(e) NA_real_)
  # Enforce column types so rows from different models always bind (list_rbind is type-strict).
  out_row <- function(method, b, se, p, n, n_out = NA_integer_, gp = NA_real_)
    tibble(model = model, method = method, n_snp = as.integer(n),
           b = as.numeric(b), se = as.numeric(se), p = as.numeric(p),
           Q_p = as.numeric(Q_p), egger_int_p = as.numeric(egg_p),
           n_outlier = as.integer(n_out), presso_global_p = as.numeric(gp))

  # MR-PRESSO reports very small p-values as strings ("<0.001"); coerce to numeric so rows bind.
  as_num <- function(x) if (is.character(x)) suppressWarnings(as.numeric(sub("^<", "", x))) else as.numeric(x)

  if (model == "MR-PRESSO") {                              # SOP: chosen when outliers seen / Q sig
    pr <- tryCatch(suppressWarnings(run_mr_presso(dat, NbDistribution = 1000,
                                                  SignifThreshold = 0.05)), error = \(e) NULL)
    if (!is.null(pr)) {
      main <- pr[[1]]$`Main MR results`
      gp   <- as_num(pr[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
      oidx <- pr[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      n_out <- if (is.null(oidx)) 0L else length(oidx)
      row <- main[main$`MR Analysis` == "Outlier-corrected", ]
      lab <- "MR-PRESSO (outlier-corrected)"
      if (nrow(row) == 0 || is.na(as_num(row$`Causal Estimate`))) {   # no outliers -> raw estimate
        row <- main[main$`MR Analysis` == "Raw", ]; lab <- "MR-PRESSO (raw; no outliers)"
      }
      return(out_row(lab, as_num(row$`Causal Estimate`), as_num(row$Sd),
                     as_num(row$`P-value`), nrow(dat), n_out, gp))
    }
    message("MR-PRESSO unavailable/failed — falling back to IVW-MRE"); model <- "IVW-MRE"
  }

  meth <- unname(tsm[model])   # single-bracket: NA for an unknown model (not an error)
  if (is.na(meth)) { message("unknown model '", model, "' — using IVW-MRE")
                     model <- "IVW-MRE"; meth <- "mr_ivw_mre" }
  res <- tryCatch(mr(dat, method_list = meth), error = \(e) NULL)
  if (is.null(res) || nrow(res) == 0) {                                # e.g. mr.raps not installed
    res <- mr(dat, method_list = "mr_ivw_mre")
    return(out_row("IVW-MRE (fallback)", res$b[1], res$se[1], res$pval[1], res$nsnp[1]))
  }
  out_row(model, res$b[1], res$se[1], res$pval[1], res$nsnp[1])
}

# One exposure x one AD outcome: pull, harmonise, run the SOP-validated model -> one estimate row.
mr_pair <- function(exp_id, exp_label, out_id, out_label, model = "IVW-MRE",
                    p = P_THRESH, r2 = R2, kb = KB) {
  inst <- get_instruments(exp_id, p, r2, kb)
  if (is.null(inst) || nrow(inst) == 0) stop("no instruments: ", exp_id)
  out <- og_retry(function() extract_outcome_data(snps = inst$SNP, outcomes = out_id))
  if (is.null(out) || nrow(out) == 0) stop("no overlap: ", out_id)
  dat <- harmonise_data(inst, out, action = 2) |> filter(mr_keep)
  if (nrow(dat) == 0) stop("nothing kept after harmonisation")
  run_naive_model(dat, model) |>
    mutate(exposure = exp_label, ad = out_label, .before = 1) |>
    rename(b_ad = b, se_ad = se, p_ad = p)
}

# Typed empty prototype so a fully-failed pull still carries the schema downstream joins need.
naive_proto <- tibble(exposure = character(), ad = character(), model = character(),
                      method = character(), n_snp = integer(), b_ad = double(),
                      se_ad = double(), p_ad = double(), Q_p = double(),
                      egger_int_p = double(), n_outlier = integer(), presso_global_p = double())
```
:::



::: {.cell}

```{.r .cell-code}
naive_grid <- expand_grid(
  exposures |> transmute(exp_id = id, exp_label = label, model = naive_model),
  tibble(out_id = unname(ad_gwas), out_label = names(ad_gwas))
)

# Persist the successful pull to CSV rather than knitr cache: an empty/failed pull is never
# written, so a later render retries once the OpenGWAS allowance is back.
naive_csv <- "exposure_ad_naive.csv"
# Re-pull if the cache is missing, predates the schema, OR was computed under different
# `naive_model` choices than the config now requests — so changing a model auto-invalidates it.
cache_ok <- local({
  if (!file.exists(naive_csv)) return(FALSE)
  hdr <- names(read_csv(naive_csv, n_max = 0, show_col_types = FALSE))
  if (!all(c("exposure", "model", "b_ad", "se_ad", "p_ad", "n_outlier", "Q_p") %in% hdr))
    return(FALSE)
  cached    <- read_csv(naive_csv, show_col_types = FALSE) |> distinct(exposure, model)
  requested <- exposures |> transmute(exposure = label, model = naive_model)
  nrow(anti_join(requested, cached, by = c("exposure", "model"))) == 0
})

naive_ad <-
  if (cache_ok) {
    read_csv(naive_csv, show_col_types = FALSE)
  } else {
    x <- pmap(naive_grid, possibly(mr_pair, otherwise = NULL)) |> compact() |> list_rbind()
    if (nrow(x) > 0) {
      write_csv(x, naive_csv); x
    } else {
      warning("Naive E→AD returned no rows — OpenGWAS is likely rate-limited or the token's ",
              "allowance is exhausted. Check ieugwasr::user(); rerun once available. The ",
              "correction below will be empty until then.")
      naive_proto
    }
  }

naive_ad |>
  select(exposure, ad, model, method, n_snp, b_ad, se_ad, p_ad, n_outlier, Q_p) |>
  mutate(across(where(is.numeric), \(x) signif(x, 3))) |>
  knitr::kable(caption = "Naive exposure → AD, per SOP-validated model (n_outlier / Q p = MR-PRESSO / heterogeneity diagnostics).")
```

::: {.cell-output-display}


Table: Naive exposure → AD, per SOP-validated model (n_outlier / Q p = MR-PRESSO / heterogeneity diagnostics).

|exposure                   |ad                                 |model     |method                        | n_snp|    b_ad|  se_ad|    p_ad| n_outlier|    Q_p|
|:--------------------------|:----------------------------------|:---------|:-----------------------------|-----:|-------:|------:|-------:|---------:|------:|
|Smoking (GSCAN, ieu-b-142) |AD — Bellenguez 2022 (proxy-incl.) |MR-PRESSO |MR-PRESSO (outlier-corrected) |    23| -0.1120| 0.0305| 0.00147|         2| 0.0000|
|Smoking (GSCAN, ieu-b-142) |AD — Kunkle 2019 IGAP (clinical)   |MR-PRESSO |MR-PRESSO (raw; no outliers)  |    21| -0.0263| 0.0708| 0.71400|         0| 0.0416|
|LDL-C                      |AD — Bellenguez 2022 (proxy-incl.) |MR-PRESSO |MR-PRESSO (outlier-corrected) |   161| -0.0718| 0.0344| 0.03870|         8| 0.0000|
|LDL-C                      |AD — Kunkle 2019 IGAP (clinical)   |MR-PRESSO |MR-PRESSO (outlier-corrected) |   149|  0.0795| 0.0654| 0.22700|        16| 0.0000|


:::
:::




Model selection (scatter/funnel inspection, MR-PRESSO vs IVW-MRE vs …) was done
when these associations were first identified; the validated choice is recorded in
`naive_model` above and applied here. This notebook only carries the correction.

## Exposure → lifespan arm (reused from Step 1)

`E → lifespan` on the same selection axis `b_SH` is fit on. Read from Step 1's
saved screen; if that file is absent, recompute with the same function.


::: {.cell}

```{.r .cell-code}
es_arm <-
  if (file.exists("exposure_survival_screen.csv")) {
    read_csv("exposure_survival_screen.csv", show_col_types = FALSE) |>
      filter(survival == SELECTION_LABEL) |>
      transmute(exposure, b_life = b, se_life = se, p_life = p)
  } else {
    warning("exposure_survival_screen.csv not found — recomputing E → lifespan (IVW-MRE).")
    life_id <- "ebi-a-GCST006697"
    exp_tbl <- exposures |> transmute(exp_id = id, exp_label = label)
    pmap(exp_tbl, \(exp_id, exp_label)
         possibly(mr_pair, NULL)(exp_id, exp_label, life_id, SELECTION_LABEL)) |>
      compact() |> list_rbind() |>
      transmute(exposure, b_life = b_ad, se_life = se_ad, p_life = p_ad)
  }

es_arm |>
  mutate(across(where(is.numeric), \(x) signif(x, 3))) |>
  knitr::kable(caption = "Exposure → lifespan (selection axis), from Step 1.")
```

::: {.cell-output-display}


Table: Exposure → lifespan (selection axis), from Step 1.

|exposure               |   b_life|  se_life|  p_life|
|:----------------------|--------:|--------:|-------:|
|chronotype.ukbb        | -0.00327| 0.043200| 0.94000|
|sleeplessness.elsworth | -0.01520| 0.051600| 0.76900|
|sleepduration.elsworth | -0.07080| 0.040800| 0.08260|
|sleepapnea.sakaue      |  0.02530| 0.021900| 0.24700|
|LDLchol.ukb.richardson |  0.06500| 0.010600| 0.00000|
|HDLchol.ukb.richardson | -0.06560| 0.008860| 0.00000|
|weight.ukb             |  0.10400| 0.010100| 0.00000|
|FBG.manning            |  0.03510| 0.017700| 0.04760|
|PBG.chen               |  0.02090| 0.012900| 0.10400|
|T2D.xue                |  0.02730| 0.004940| 0.00000|
|BMI.elsworth           |  0.13600| 0.011100| 0.00000|
|bloodglucose.barton    |  0.01430| 0.010700| 0.18100|
|HbA1C.soranzo          | -0.01800| 0.037600| 0.63200|
|schizophrenia.pgc      |  0.00582| 0.005380| 0.28000|
|MDD.pgc                | -0.00869| 0.023000| 0.70600|
|bipolar.pgc            | -0.02210| 0.007080| 0.00177|
|menopause.ukb          | -0.00796| 0.084700| 0.92500|
|SBP.tang               |  0.01000| 0.000812| 0.00000|
|DBP.tang               |  0.01640| 0.001410| 0.00000|
|EA.yearsofed           | -0.03480| 0.003020| 0.00000|
|EA.collegecomplete     | -0.38400| 0.032500| 0.00000|
|CVD.loh                |  0.60300| 0.030100| 0.00000|
|cig.liu                |  0.10200| 0.017700| 0.00000|
|alcohol.liu            |  0.08870| 0.039700| 0.02560|
|fatmass.elsworth       |  0.12700| 0.011800| 0.00000|


:::
:::


## Selection slope b_SH

One genome-wide, APOE-excluded SlopeHunter slope per AD GWAS. Bellenguez × this
lifespan is reused from the smoking pipeline; a missing fit (e.g. Kunkle before
you run the driver) yields `NA` and its exposures are reported as a gap rather
than silently left uncorrected.


::: {.cell}

```{.r .cell-code}
load_bsh <- function(ad_id) {
  src <- bsh_sources |> filter(ad_id == !!ad_id)
  if (nrow(src) == 0 || !file.exists(src$fits_csv[1]))
    return(tibble(ad_id = ad_id, b_SH = NA_real_, se_SH = NA_real_, n_fit = NA_integer_))
  read_csv(src$fits_csv[1], show_col_types = FALSE) |>
    filter(set == src$fit_set[1]) |>
    transmute(ad_id = ad_id, b_SH, se_SH, n_fit = dplyr::coalesce(n_fit, NA_integer_)) |>
    slice(1)
}

bsh <- map(unname(ad_gwas), load_bsh) |> list_rbind() |>
  left_join(tibble(ad_id = unname(ad_gwas), ad = names(ad_gwas)), by = "ad_id")

bsh |>
  transmute(`AD GWAS` = ad, b_SH = signif(b_SH, 3), se_SH = signif(se_SH, 3), n_fit) |>
  knitr::kable(caption = "SlopeHunter selection slope b_SH (genome-wide, APOE excluded).")
```

::: {.cell-output-display}


Table: SlopeHunter selection slope b_SH (genome-wide, APOE excluded).

|AD GWAS                            |   b_SH|  se_SH| n_fit|
|:----------------------------------|------:|------:|-----:|
|AD — Bellenguez 2022 (proxy-incl.) | -0.944| 0.0701|  2250|
|AD — Kunkle 2019 IGAP (clinical)   |  1.100| 0.1120|  2250|


:::
:::


To fit `b_SH` for an AD GWAS that has no cached fit yet (needs full sumstats +
plink2 clumping — see `smoking/scripts/fetch_sumstats.sh`):


::: {.cell}

```{.r .cell-code}
# Turnkey: `bash scripts/fetch_kunkle_sumstats.sh` downloads + normalizes the Kunkle and
# lifespan sumstats and the LD panel, then prints the prep + `plink2 --clump` commands. After
# clumping, this single call preps, filters to the clumped set, fits, and writes the fits CSV:
dir.create("R/fits", showWarnings = FALSE, recursive = TRUE)
build_selection_slope(
  ad_file      = "data/cache/sumstats/AD_kunkle_harmonised.tsv.gz",
  life_file    = "data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz",
  clumped_snps = "data/cache/sumstats/kunkle_clumped.clumps",  # plink2 --clump output
  out_csv      = "R/fits/kunkle_lifespan_slopehunter_fits.csv")
```
:::


::: callout-warning
**Verify the sign of a new `b_SH` before trusting it.** The current Kunkle fit is
`b_SH = +1.10`, opposite to the validated Bellenguez `−0.944`. Because AD → shorter
lifespan holds regardless of ascertainment, an opposite sign is a red flag for an
allele-orientation artefact in the AD sumstats — and it flips the direction of every
Kunkle correction. Confirm before reporting the Kunkle column:
:::


::: {.cell}

```{.r .cell-code}
# Aligned AD betas for the same disease must correlate POSITIVELY across the two AD GWAS.
# A negative correlation => the test GWAS's effect alleles are inverted => its b_SH is artefactual.
m_bell <- prep_selection_merge("smoking/data/cache/sumstats/AD_bellenguez_harmonised.tsv.gz",
                               "data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz")
m_kunk <- prep_selection_merge("data/cache/sumstats/AD_kunkle_harmonised.tsv.gz",
                               "data/cache/sumstats/LIFESPAN_pilling_harmonised.tsv.gz")
ad_orientation_check(m_bell, m_kunk, "Bellenguez", "Kunkle")   # verdict: consistent vs FLIPPED
```
:::


## Corrected exposure → AD


::: {.cell}

```{.r .cell-code}
joined <- naive_ad |>
  left_join(es_arm, by = "exposure") |>
  left_join(bsh |> select(ad, b_SH, se_SH), by = "ad")

# correct_estimate() is vectorized, so call it once on whole columns. This is NA-safe (b_SH NA
# -> b_corr NA) and empty-safe (0 rows in -> a 0-row tibble that still carries b_corr etc.),
# unlike a per-row list-column + unnest() which cannot materialize columns from zero elements.
corrected <- bind_cols(
  joined,
  correct_estimate(joined$b_ad, joined$se_ad, joined$b_life, joined$se_life,
                   joined$b_SH, joined$se_SH)
) |>
  mutate(
    attenuation = 1 - abs(b_corr) / abs(b_uncorr),
    verdict = case_when(
      is.na(b_corr)                                   ~ "no b_SH — run the fit for this AD GWAS",
      sign(b_corr) != sign(b_uncorr) & p_corr < 0.05  ~ "REVERSED after correction",
      p_ad >= 0.05 & p_corr < 0.05                    ~ "UNMASKED — null naive, significant after correction",
      p_ad < 0.05 & p_corr >= 0.05                    ~ "ABOLISHED — signal explained by survival bias",
      p_corr < 0.05                                   ~ "ROBUST — survives correction",
      TRUE                                            ~ "inconclusive"
    )
  )
```
:::



::: {.cell}

```{.r .cell-code}
corrected |>
  transmute(
    exposure, `AD GWAS` = ad, model, n_snp,
    naive = signif(b_ad, 3), `p naive` = signif(p_ad, 2),
    b_SH = signif(b_SH, 3), `E→life` = signif(b_life, 3),
    corrected = signif(b_corr, 3),
    `corr 95% CI` = ifelse(is.na(b_corr), NA,
                           sprintf("[%.3f, %.3f]", ci_lo, ci_hi)),
    `p corr` = signif(p_corr, 2),
    `% attenuated` = ifelse(is.na(attenuation), NA, round(100 * attenuation)),
    verdict
  ) |>
  knitr::kable(caption = paste0(
    "Naive vs SlopeHunter-corrected exposure → AD. ",
    "corrected = naive − b_SH × (E→lifespan)."))
```

::: {.cell-output-display}


Table: Naive vs SlopeHunter-corrected exposure → AD. corrected = naive − b_SH × (E→lifespan).

|exposure                   |AD GWAS                            |model     | n_snp|   naive| p naive|   b_SH| E→life| corrected|corr 95% CI | p corr|% attenuated |verdict                                |
|:--------------------------|:----------------------------------|:---------|-----:|-------:|-------:|------:|------:|---------:|:-----------|------:|:------------|:--------------------------------------|
|Smoking (GSCAN, ieu-b-142) |AD — Bellenguez 2022 (proxy-incl.) |MR-PRESSO |    23| -0.1120|  0.0015| -0.944|     NA|        NA|NA          |     NA|NA           |no b_SH — run the fit for this AD GWAS |
|Smoking (GSCAN, ieu-b-142) |AD — Kunkle 2019 IGAP (clinical)   |MR-PRESSO |    21| -0.0263|  0.7100|  1.100|     NA|        NA|NA          |     NA|NA           |no b_SH — run the fit for this AD GWAS |
|LDL-C                      |AD — Bellenguez 2022 (proxy-incl.) |MR-PRESSO |   161| -0.0718|  0.0390| -0.944|     NA|        NA|NA          |     NA|NA           |no b_SH — run the fit for this AD GWAS |
|LDL-C                      |AD — Kunkle 2019 IGAP (clinical)   |MR-PRESSO |   149|  0.0795|  0.2300|  1.100|     NA|        NA|NA          |     NA|NA           |no b_SH — run the fit for this AD GWAS |


:::
:::


## Naive vs corrected — forest


::: {.cell}

```{.r .cell-code}
plot_df <- corrected |> filter(!is.na(b_corr))

if (nrow(plot_df) == 0) {
  message("No corrected estimates to plot yet (naive pull empty or b_SH unavailable).")
} else {
  plot_df |>
    mutate(lab = paste(exposure, ad, sep = " → ")) |>
    select(lab, Naive = b_ad, Corrected = b_corr,
           se_n = se_ad, se_c = se_corr) |>
    pivot_longer(c(Naive, Corrected), names_to = "estimate", values_to = "b") |>
    mutate(se = if_else(estimate == "Naive", se_n, se_c)) |>
    ggplot(aes(b, lab, colour = estimate)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                    position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c(Naive = "#ffcb05", Corrected = "#00274c")) +
    labs(x = "effect on AD (log-OR)", y = NULL, colour = NULL,
         title = "Survival-collider correction of exposure → AD") +
    theme_minimal(base_size = 12)
}
```
:::


## How to read this

- **ABOLISHED** — the naive effect is explained by the survival collider: after
  removing `b_SH · (E→lifespan)` the corrected CI covers 0. Do not treat the
  naive `E → AD` effect as causal/druggable (this was the smoking→AD verdict).
- **ROBUST** — significant before *and* after; the collider does not account for
  it (though see caveats).
- **UNMASKED** — the naive was null but the corrected effect is significant. The
  collider was *suppressing* a real effect: `naive = true + bias`, and here the
  bias ≈ −true, so the two cancelled to ~0 until correction removed the bias.
  A null naive never means "no effect" — the collider can hide one.
- **REVERSED** — correction flips the sign; the induced bias dominated the naive.

**Trust caveat for UNMASKED/REVERSED.** These "revealed" effects are only as
reliable as `b_SH` and `E→lifespan`. A wrong selection slope (e.g. a mis-oriented
`b_SH`) can *manufacture* a significant corrected effect out of a true null just
as easily as it can reveal a real one — so verify the inputs (allele orientation,
`b_SH` sign) before believing an UNMASKED result.

The `induced_bias` direction from Step 1 predicts *which way* the naive estimate
is pulled; this step quantifies *how much*. They should agree in sign: a
"spurious PROTECTION" flag in Step 1 means correction should move the naive
`E → AD` estimate in the harmful (less protective) direction.

Caveats carried from Step 1 and the smoking analysis:

1. **b_SH is per AD GWAS.** Bellenguez (proxy-inclusive) and Kunkle (clinical)
   can carry different selection slopes; fit each rather than reusing one.
2. **The correction cannot separate collider-induced association from a genuine
   effect collinear with the mortality axis** — strong evidence, not proof.
3. **`E → lifespan` uses kin-proxy parental lifespan** (~halved per-allele
   effects), so the subtracted term is conservative.


::: {.cell}

```{.r .cell-code}
write_csv(corrected, "exposure_ad_indexevent_corrected.csv")
```
:::

