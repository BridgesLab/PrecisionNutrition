---
title: "MR Analyses of Secondary Outcomes for Total Cholesterol on Calcium Homeostasis"
author: "Dave Bridges"
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
execute:
  echo: true
  warning: false
bibliography: references.bib
---


::: {.cell}

:::


## Purpose

To test if SNPs for total cholesterol GWAS identified using UK Biobank relate
to other mechanistic or pathological outcomes related to calcium homeostasis
and bone health. This script can be found in /Users/davebrid/Documents/GitHub/PrecisionNutrition/Human Genetics and was most recently
run on Thu Apr 30 08:06:05 2026.

This is a revised version of an earlier analysis. The previous version used
local PheWeb summary statistics (GEFOS 2012/2015) for bone outcomes and
encountered substantial instrument loss during harmonisation — only 84 of 370
total cholesterol instruments survived for lumbar spine BMD, with no usable
overlap for the 2012 femoral neck GWAS. This version replaces those local
files with OpenGWAS API queries against larger, more recent BMD and fracture
GWAS, using LD proxies (r²>0.8) to recover instruments that aren't directly
present in each outcome dataset.

## Data Entry


::: {.cell}

```{.r .cell-code}
instruments.tc.file <- 'Total Cholesterol Instruments from UKBB.csv'

instruments.tc <- read_csv(instruments.tc.file, show_col_types = FALSE) |>
  dplyr::rename(
    SNP                       = SP2,
    beta.exposure             = BETA,
    se.exposure               = SE,
    effect_allele.exposure    = EA,
    other_allele.exposure     = OA,
    pval.exposure             = P,
    eaf.exposure              = ALT_FREQS,
    samplesize.exposure       = N_exposure
  ) |>
  mutate(id.exposure = "Total Cholesterol (UK Biobank)",
         exposure    = "Total Cholesterol (UK Biobank)")
```
:::


We used 370 SNPs as instruments for total cholesterol from
UK Biobank. These are found in the Total Cholesterol Instruments from UKBB.csv datafile.

### Converting Instrument Identifiers to rsIDs

The instrument file uses chromosome-position-allele IDs (e.g. `1:55505647:C:T`)
but OpenGWAS queries require rsIDs. We use `SNPlocs.Hsapiens.dbSNP144.GRCh37`
to look up rsIDs by genomic position.


::: {.cell}

```{.r .cell-code}
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(GenomicRanges)

# Force tidyverse versions back to top of search path after Bioconductor masking
rename  <- dplyr::rename
first   <- dplyr::first
compact <- purrr::compact

snpdb <- SNPlocs.Hsapiens.dbSNP144.GRCh37

tc.positions <- instruments.tc |>
  separate(SNP, into = c("CHR", "BP", "ALLELE0", "ALLELE1"),
           sep = ":", remove = FALSE) |>
  mutate(CHR = as.integer(CHR), BP = as.integer(BP))

get_rsids_by_chr <- function(chr, positions, snpdb) {
  pos_ranges <- GPos(seqnames = chr, pos = positions)
  snps <- snpsByOverlaps(snpdb, pos_ranges)
  data.frame(
    CHR  = as.integer(chr),
    BP   = as.integer(pos(snps)),
    RSID = snps$RefSNP_id,
    stringsAsFactors = FALSE
  )
}

rsid_list <- list()
for (chr in unique(tc.positions$CHR)) {
  chr_data <- tc.positions |> filter(CHR == chr)
  rsid_list[[chr]] <- get_rsids_by_chr(chr, chr_data$BP, snpdb)
  message("Processed chr", chr)
}
all_rsids <- bind_rows(rsid_list) %>%
  mutate(CHR = as.integer(CHR), BP = as.integer(BP))

# Build rsID-keyed instruments in one pipe to avoid the dplyr/Bioc masking bug
instruments.tc.rsid <- instruments.tc %>%
  separate(SNP, into = c("CHR", "BP", "ALLELE0", "ALLELE1"),
           sep = ":", remove = TRUE) %>%
  mutate(CHR = as.integer(CHR), BP = as.integer(BP)) %>%
  left_join(all_rsids, by = c("CHR", "BP")) %>%
  filter(!is.na(RSID)) %>%
  dplyr::select(
    effect_allele.exposure, other_allele.exposure,
    beta.exposure, se.exposure, pval.exposure, eaf.exposure,
    samplesize.exposure,
    id.exposure, exposure,
    RSID
  ) %>%
  dplyr::rename(SNP = RSID)

cat("Original instruments:              ", nrow(instruments.tc), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Original instruments:               370 
```


:::

```{.r .cell-code}
cat("After rsID lookup (excluding NAs): ", nrow(instruments.tc.rsid), "\n")
```

::: {.cell-output .cell-output-stdout}

```
After rsID lookup (excluding NAs):  360 
```


:::
:::


## Mechanistic Outcomes

### Vitamin D Levels

This analysis tests the hypothesis that total cholesterol impacts
25-hydroxyvitamin D levels, as cholesterol is a precursor for vitamin D
synthesis in the skin. This could positively impact calcium levels indirectly
via increased vitamin D. This analysis uses MGI-BioVU LabWAS data
[@goldsteinLabWASNovelFindings2020]. The vitamin D GWAS is retained from the
original local file because it is a non-routine biomarker with limited
coverage on OpenGWAS.


::: {.cell}

```{.r .cell-code}
gwas.vitd.file <- 'PheWeb Summary Statistics/phenocode-Vit-D.tsv.gz'
samplesize.outcome.vitd <- 12250

gwas.vitd <- read_tsv(gwas.vitd.file, show_col_types = FALSE) |>
  mutate(ID = paste(chrom, pos, ref, alt, sep = ":")) |>
  dplyr::rename(
    SNP                   = ID,
    beta.outcome          = beta,
    se.outcome            = sebeta,
    effect_allele.outcome = alt,
    other_allele.outcome  = ref,
    pval.outcome          = pval,
    eaf.outcome           = maf,
  ) |>
  mutate(
    id.outcome         = "Vitamin D (MGI-BioVU LabWAS)",
    outcome            = "Vitamin D (MGI-BioVU LabWAS)",
    samplesize.outcome = samplesize.outcome.vitd
  )

library(TwoSampleMR)

# Vitamin D analysis uses the original chr:pos:a:b SNP IDs (not rsIDs)
vitd.data <- harmonise_data(instruments.tc, gwas.vitd, action = 2)
vitd.data_steiger <- steiger_filtering(vitd.data)

# Pre-harmonization instrument metrics (all 370 input SNPs)
pre_harm_metrics_vitd <- instruments.tc %>%
  mutate(
    R2.exposure = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
    F.exposure  = (R2.exposure * (samplesize.exposure - 2)) / (1 - R2.exposure)
  )

pre_harm_summary_vitd <- pre_harm_metrics_vitd %>%
  summarise(
    num_snps            = n(),
    samplesize.exposure = dplyr::first(samplesize.exposure),
    cumulative_R2       = sum(R2.exposure, na.rm = TRUE),
    mean_F              = mean(F.exposure, na.rm = TRUE),
    median_F            = median(F.exposure, na.rm = TRUE),
    mean_maf            = mean(eaf.exposure, na.rm = TRUE),
    mean_beta           = mean(abs(beta.exposure), na.rm = TRUE)
  ) |>
  mutate(
    overall_F = (cumulative_R2 * (samplesize.exposure - num_snps - 1)) /
                ((1 - cumulative_R2) * num_snps)
  )

# Post-harmonization instrument metrics (SNPs that survived harmonization+Steiger)
vitd.data.annot <- vitd.data_steiger %>%
  mutate(
    R2.exposure = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
    F.exposure  = (R2.exposure * (samplesize.exposure - 2)) / (1 - R2.exposure)
  )

vitd.exposure.summary <- vitd.data.annot %>%
  summarise(
    num_snps            = n(),
    samplesize.exposure = dplyr::first(samplesize.exposure),
    cumulative_R2       = sum(R2.exposure, na.rm = TRUE),
    mean_F              = mean(F.exposure, na.rm = TRUE),
    median_F            = median(F.exposure, na.rm = TRUE),
    mean_maf            = mean(eaf.exposure, na.rm = TRUE),
    mean_beta           = mean(abs(beta.exposure), na.rm = TRUE)
  ) |>
  mutate(
    overall_F = (cumulative_R2 * (samplesize.exposure - num_snps - 1)) /
                ((1 - cumulative_R2) * num_snps)
  )

# Write instrument files
pre_harm_summary_vitd %>%
  write_csv("Instrument Metrics - Total Cholesterol for Vitamin D - Pre-Harmonization.csv")
vitd.exposure.summary %>%
  write_csv("Instrument Metrics - Total Cholesterol for Vitamin D - Post-Harmonization.csv")
vitd.data.annot %>%
  write_csv("Total Cholesterol Instruments for Vitamin D.csv")

bind_rows(
  pre_harm_summary_vitd  %>% mutate(Stage = "Pre-Harmonization"),
  vitd.exposure.summary  %>% mutate(Stage = "Post-Harmonization")
) %>%
  dplyr::select(Stage, everything()) %>%
  kable(caption = "Total cholesterol instruments before and after harmonisation for Vitamin D analysis",
        digits = c(NA, 0, 0, 4, 1, 1, 4, 4, 1))
```

::: {.cell-output-display}


Table: Total cholesterol instruments before and after harmonisation for Vitamin D analysis

|Stage              | num_snps| samplesize.exposure| cumulative_R2| mean_F| median_F| mean_maf| mean_beta| overall_F|
|:------------------|--------:|-------------------:|-------------:|------:|--------:|--------:|---------:|---------:|
|Pre-Harmonization  |      370|              420607|        0.1061|  120.8|     48.9|   0.3171|    0.0309|     134.8|
|Post-Harmonization |      285|              420607|        0.0972|  143.8|     53.5|   0.3498|    0.0301|     158.8|


:::

```{.r .cell-code}
vitd.mr <- mr(vitd.data_steiger,
              method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_raps",
                              "mr_egger_regression",
                              "mr_weighted_median", "mr_weighted_mode"))

vitd.mr |>
  dplyr::select(-starts_with('id')) |>
  kable(caption = "MR Results for Total Cholesterol - Vitamin D Analysis",
        digits = c(0, 0, 0, 0, 3, 3, 99))
```

::: {.cell-output-display}


Table: MR Results for Total Cholesterol - Vitamin D Analysis

|outcome                      |exposure                       |method                                                    | nsnp|      b|    se|       pval|
|:----------------------------|:------------------------------|:---------------------------------------------------------|----:|------:|-----:|----------:|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Inverse variance weighted (multiplicative random effects) |  280| -0.063| 0.033| 0.05362791|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Inverse variance weighted (fixed effects)                 |  280| -0.063| 0.029| 0.02795864|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Robust adjusted profile score (RAPS)                      |  280| -0.067| 0.034| 0.04513899|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |MR Egger                                                  |  280| -0.070| 0.053| 0.19017784|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted median                                           |  280| -0.005| 0.051| 0.92114216|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted mode                                             |  280| -0.012| 0.055| 0.83388196|


:::

```{.r .cell-code}
vitd.pleio <- mr_pleiotropy_test(vitd.data_steiger)
vitd.pleio |>
  dplyr::select(-starts_with('id')) |>
  kable(caption = "MR Pleiotropy Results for Total Cholesterol - Vitamin D Analysis")
```

::: {.cell-output-display}


Table: MR Pleiotropy Results for Total Cholesterol - Vitamin D Analysis

|outcome                      |exposure                       | egger_intercept|        se|      pval|
|:----------------------------|:------------------------------|---------------:|---------:|---------:|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |       0.0002851| 0.0017118| 0.8678594|


:::

```{.r .cell-code}
vitd.het <- mr_heterogeneity(vitd.data_steiger) |>
  mutate(I2 = pmax(0, (Q - Q_df) / Q) * 100)
vitd.het |>
  dplyr::select(-starts_with('id')) |>
  kable(caption = "MR Heterogeneity Results for Total Cholesterol - Vitamin D Analysis",
        digits = c(0, 0, 0, 3, 3, 99, 1))
```

::: {.cell-output-display}


Table: MR Heterogeneity Results for Total Cholesterol - Vitamin D Analysis

|outcome                      |exposure                       |method                    |       Q| Q_df|       Q_pval|   I2|
|:----------------------------|:------------------------------|:-------------------------|-------:|----:|------------:|----:|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |MR Egger                  | 361.847|  278| 0.0005223738| 23.2|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Inverse variance weighted | 361.883|  279| 0.0005999283| 22.9|


:::
:::


### Vitamin D — Diagnostic Plots

#### Scatter Plot


::: {.cell}

```{.r .cell-code}
ggplot(vitd.data_steiger, aes(x = beta.exposure, y = beta.outcome)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome,
                    ymax = beta.outcome + 1.96 * se.outcome),
                alpha = 0.5) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96 * se.exposure,
                    xmax = beta.exposure + 1.96 * se.exposure),
                alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 16) +
  labs(x = "Exposure Estimate (Total Cholesterol)",
       y = "Outcome Estimate (Vitamin D)",
       title = "")
```

::: {.cell-output-display}
![](figures/vitd-scatter-1.png){width=672}
:::
:::


#### Funnel Plot


::: {.cell}

```{.r .cell-code}
vitd.single_snp <- mr_singlesnp(vitd.data_steiger)

vitd.ivw_beta <- vitd.mr |>
  filter(method == "Inverse variance weighted (multiplicative random effects)") |>
  pull(b)

vitd.y_max <- max(vitd.single_snp$se^{-1}, na.rm = TRUE) * 1.1
vitd.precision_grid <- seq(0, vitd.y_max, length.out = 1000)
vitd.bounds_df <- data.frame(
  precision = vitd.precision_grid,
  lower     = vitd.ivw_beta - 1.96 / vitd.precision_grid,
  upper     = vitd.ivw_beta + 1.96 / vitd.precision_grid
)

ggplot(vitd.single_snp, aes(x = b, y = 1/se)) +
  geom_point(size = 1) +
  geom_vline(xintercept = vitd.ivw_beta, linetype = "solid",
             color = "#ff7f0e", size = 1) +
  geom_line(data = vitd.bounds_df,
            aes(x = lower, y = precision), linetype = "dashed") +
  geom_line(data = vitd.bounds_df,
            aes(x = upper, y = precision), linetype = "dashed") +
  labs(x = "Estimate (Beta-IVW)", y = "Precision (1/Standard Error)",
       title = "") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, vitd.y_max),
                  xlim = c(min(vitd.single_snp$b, na.rm = TRUE),
                           max(vitd.single_snp$b, na.rm = TRUE)))
```

::: {.cell-output-display}
![](figures/vitd-funnel-1.png){width=672}
:::
:::


#### Leave-One-Out Analysis


::: {.cell}

```{.r .cell-code}
vitd.loo_res <- mr_leaveoneout(vitd.data_steiger)

vitd.loo_res |>
  mutate(diff = b - filter(vitd.mr,
                           method == "Inverse variance weighted (multiplicative random effects)")$b) |>
  arrange(-abs(diff)) |>
  head() |>
  dplyr::select(SNP, diff, b, se, p) |>
  kable(caption = "Leave-One-Out Results for Vitamin D Analysis (IVW method) for influential SNPs",
        digits = c(0, 5, 5, 5, 5))
```

::: {.cell-output-display}


Table: Leave-One-Out Results for Vitamin D Analysis (IVW method) for influential SNPs

|SNP             |     diff|        b|      se|       p|
|:---------------|--------:|--------:|-------:|-------:|
|1:63112320:C:G  | -0.01110| -0.07410| 0.03269| 0.02342|
|9:136154168:T:C |  0.00879| -0.05421| 0.03258| 0.09612|
|11:61569306:C:G | -0.00797| -0.07097| 0.03249| 0.02894|
|19:45349369:T:C | -0.00789| -0.07089| 0.03320| 0.03275|
|15:58680954:T:C |  0.00690| -0.05610| 0.03266| 0.08590|
|15:58723426:A:G |  0.00669| -0.05631| 0.03291| 0.08709|


:::

```{.r .cell-code}
ggplot(vitd.loo_res, aes(x = reorder(SNP, -b), y = b)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se),
                width = 0.01, alpha = 0.5) +
  coord_flip() +
  labs(x = "SNP Removed", y = "Estimate (Beta-IVW; leave-one-out)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_text(size = 1))
```

::: {.cell-output-display}
![](figures/vitd-loo-1.png){width=672}
:::
:::


Leave-one-out analyses suggested that no individual SNP had a relatively large
influence on the IVW estimate, supporting the robustness of the causal
inference for vitamin D.

#### MR-PRESSO Analysis

MR-PRESSO (Pleiotropy RESidual Sum and Outlier) tests whether the IVW estimate
is driven by a small number of pleiotropic outliers and provides an
outlier-corrected estimate.


::: {.cell}

```{.r .cell-code}
library(MRPRESSO)

set.seed(2026)
vitd.mrpresso <- tryCatch(
  mr_presso(BetaOutcome     = "beta.outcome",
            BetaExposure    = "beta.exposure",
            SdOutcome       = "se.outcome",
            SdExposure      = "se.exposure",
            OUTLIERtest     = TRUE,
            DISTORTIONtest  = TRUE,
            data            = vitd.data_steiger,
            NbDistribution  = 50000,
            SignifThreshold = 0.05),
  error = function(e) { message("MR-PRESSO failed: ", e$message); NULL }
)

if (!is.null(vitd.mrpresso)) {
  vitd.mrpresso$`Main MR results` |>
    kable(caption = "MR-PRESSO results for Vitamin D Analysis",
          digits = c(0, 4, 4, 4, 0, 99))

  cat("Global heterogeneity test p-value:",
      vitd.mrpresso$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
}
```

::: {.cell-output .cell-output-stdout}

```
Global heterogeneity test p-value: 3e-04 
```


:::

```{.r .cell-code}
# Format MR-PRESSO results as MR-style rows for the combined summary
extract_mrpresso_rows <- function(mrpresso_obj, exposure_label, outcome_label, n_snp) {
  if (is.null(mrpresso_obj)) return(NULL)
  res <- mrpresso_obj$`Main MR results`
  global_p_raw <- mrpresso_obj$`MR-PRESSO results`$`Global Test`$Pvalue
  
  # Coerce to numeric — handle "<0.001" style strings
  global_p <- if (is.character(global_p_raw)) {
    as.numeric(str_remove(global_p_raw, "<"))
  } else {
    as.numeric(global_p_raw)
  }

  rows <- res %>%
    mutate(
      method   = if_else(`MR Analysis` == "Raw",
                         "MR-PRESSO (Raw)", "MR-PRESSO (Corrected)"),
      exposure = exposure_label,
      outcome  = outcome_label,
      nsnp     = as.integer(n_snp),
      b        = as.numeric(`Causal Estimate`),
      se       = as.numeric(Sd),
      pval     = as.numeric(`P-value`)
    ) %>%
    dplyr::select(exposure, outcome, method, nsnp, b, se, pval)

  rows$global_pval <- global_p
  rows
}

vitd.mrpresso_rows <- extract_mrpresso_rows(
  vitd.mrpresso,
  "Total Cholesterol (UK Biobank)",
  "Vitamin D (MGI-BioVU LabWAS)",
  nrow(vitd.data_steiger)
)
```
:::



## Pathological Outcomes

### Bone Outcome GWAS Selection

The previous analysis used local PheWeb summary statistics from older GEFOS
releases (Estrada 2012; Zheng 2015). These files had three problems: small
sample sizes (32k–53k participants), build-37 marker IDs requiring manual
rsID conversion, and dramatic instrument loss during harmonisation. To
address these issues, this version queries three more recent BMD and fracture
GWAS directly via the OpenGWAS API, using LD proxies (r²>0.8) to recover
SNPs that are not directly present in each outcome dataset.

The three bone-related outcomes are treated as primary analyses because
they ask different but complementary questions:

- **Heel BMD (Morris 2019, UKB, n=426,824)**: highest-powered BMD GWAS
  available, derived from quantitative ultrasound; cortical-bone-dominated
- **Femoral neck BMD (Zheng 2015, GEFOS, n=32,735)**: DXA-based, lower
  powered legacy comparator
- **Fractures (Dönertaş 2021, UKB, n=484,598)**: clinical endpoint integrating
  bone density, geometry, and fall propensity

### Querying OpenGWAS for BMD and Fracture Outcomes


::: {.cell}

```{.r .cell-code}
library(ieugwasr)

bone_outcome_gwas <- tribble(
  ~label,                              ~gwas_id,              ~n,       ~year, ~pmid,
  "Heel BMD (Morris 2019, UKB)",        "ebi-a-GCST006979",    426824,   2019,  30598549,
  "Femoral neck BMD (Zheng 2015)",      "ieu-a-980",           32735,    2015,  26367794,
  "Fractures (Dönertaş 2021, UKB)",     "ebi-a-GCST90038703",  484598,   2021,  34187969
)

kable(bone_outcome_gwas,
      caption = "Bone outcome GWAS queried via OpenGWAS")
```

::: {.cell-output-display}


Table: Bone outcome GWAS queried via OpenGWAS

|label                          |gwas_id            |      n| year|     pmid|
|:------------------------------|:------------------|------:|----:|--------:|
|Heel BMD (Morris 2019, UKB)    |ebi-a-GCST006979   | 426824| 2019| 30598549|
|Femoral neck BMD (Zheng 2015)  |ieu-a-980          |  32735| 2015| 26367794|
|Fractures (Dönertaş 2021, UKB) |ebi-a-GCST90038703 | 484598| 2021| 34187969|


:::

```{.r .cell-code}
fetch_bone_outcome <- function(gwas_id, label) {
  cat("  Fetching", label, "(", gwas_id, ")...\n")

  outcome_dat <- tryCatch(
    extract_outcome_data(
      snps          = instruments.tc.rsid$SNP,
      outcomes      = gwas_id,
      proxies       = TRUE,
      rsq           = 0.8,
      align_alleles = 1,
      palindromes   = 1,
      maf_threshold = 0.01
    ),
    error = function(e) { cat("    Failed:", e$message, "\n"); NULL }
  )

  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) return(NULL)

  outcome_dat %>%
    mutate(outcome = label, id.outcome = label)
}

cat("Querying OpenGWAS for", nrow(bone_outcome_gwas), "bone outcomes...\n")
```

::: {.cell-output .cell-output-stdout}

```
Querying OpenGWAS for 3 bone outcomes...
```


:::

```{.r .cell-code}
bone_outcomes <- map2(bone_outcome_gwas$gwas_id,
                      bone_outcome_gwas$label,
                      fetch_bone_outcome) %>%
  set_names(bone_outcome_gwas$label) %>%
  purrr::compact()
```

::: {.cell-output .cell-output-stdout}

```
  Fetching Heel BMD (Morris 2019, UKB) ( ebi-a-GCST006979 )...
```


:::

::: {.cell-output .cell-output-stdout}

```
  Fetching Femoral neck BMD (Zheng 2015) ( ieu-a-980 )...
```


:::

::: {.cell-output .cell-output-stdout}

```
  Fetching Fractures (Dönertaş 2021, UKB) ( ebi-a-GCST90038703 )...
```


:::

```{.r .cell-code}
map_dfr(bone_outcomes, ~ tibble(
  n_snps_returned = nrow(.x),
  n_via_proxy     = sum(.x$proxy.outcome == TRUE, na.rm = TRUE)
), .id = "outcome") |>
  mutate(
    n_instruments = nrow(instruments.tc.rsid),
    pct_recovered = round(100 * n_snps_returned / n_instruments, 1)
  ) |>
  kable(caption = "Outcome SNP recovery via OpenGWAS")
```

::: {.cell-output-display}


Table: Outcome SNP recovery via OpenGWAS

|outcome                        | n_snps_returned| n_via_proxy| n_instruments| pct_recovered|
|:------------------------------|---------------:|-----------:|-------------:|-------------:|
|Heel BMD (Morris 2019, UKB)    |             348|          11|           360|          96.7|
|Femoral neck BMD (Zheng 2015)  |             350|          88|           360|          97.2|
|Fractures (Dönertaş 2021, UKB) |             356|           0|           360|          98.9|


:::
:::


### Running MR Across Bone Outcomes


::: {.cell}

```{.r .cell-code}
instruments.tc.for.bone <- instruments.tc.rsid

# Pre-harmonization summary for the rsID-keyed instrument set used by all bone analyses
pre_harm_metrics_bone <- instruments.tc.rsid %>%
  mutate(
    R2.exposure = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
    F.exposure  = (R2.exposure * (samplesize.exposure - 2)) / (1 - R2.exposure)
  )

pre_harm_summary_bone <- pre_harm_metrics_bone %>%
  summarise(
    num_snps            = n(),
    samplesize.exposure = dplyr::first(samplesize.exposure),
    cumulative_R2       = sum(R2.exposure, na.rm = TRUE),
    mean_F              = mean(F.exposure, na.rm = TRUE),
    median_F            = median(F.exposure, na.rm = TRUE),
    mean_maf            = mean(eaf.exposure, na.rm = TRUE),
    mean_beta           = mean(abs(beta.exposure), na.rm = TRUE)
  ) |>
  mutate(
    overall_F = (cumulative_R2 * (samplesize.exposure - num_snps - 1)) /
                ((1 - cumulative_R2) * num_snps)
  )

# Map outcome labels -> safe filenames
safe_label <- function(label) {
  label %>%
    str_replace_all("[^A-Za-z0-9 ]", "") %>%
    str_squish()
}

run_bone_mr <- function(outcome_dat, label) {

  harm <- harmonise_data(instruments.tc.for.bone, outcome_dat, action = 2)
  harm_steiger <- steiger_filtering(harm)

  if (nrow(harm_steiger) == 0) {
    cat("  ", label, ": no SNPs retained after harmonisation\n")
    return(NULL)
  }

  inst_summary <- harm_steiger %>%
    mutate(R2 = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
           F  = (R2 * (samplesize.exposure - 2)) / (1 - R2)) %>%
    summarise(
      outcome             = label,
      num_snps            = n(),
      samplesize.exposure = dplyr::first(samplesize.exposure),
      cumulative_R2       = sum(R2, na.rm = TRUE),
      mean_F              = mean(F, na.rm = TRUE),
      median_F            = median(F, na.rm = TRUE),
      mean_maf            = mean(eaf.exposure, na.rm = TRUE),
      mean_beta           = mean(abs(beta.exposure), na.rm = TRUE)
    ) %>%
    mutate(
      overall_F = (cumulative_R2 * (samplesize.exposure - num_snps - 1)) /
                  ((1 - cumulative_R2) * num_snps)
    )

  mr_res <- mr(harm_steiger,
               method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_raps",
                               "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode")) %>%
    mutate(outcome = label)

  pleio <- tryCatch(
    mr_pleiotropy_test(harm_steiger) %>% mutate(outcome = label),
    error = function(e) NULL
  )

  het <- tryCatch(
    mr_heterogeneity(harm_steiger) %>%
      mutate(outcome = label, I2 = pmax(0, (Q - Q_df) / Q) * 100),
    error = function(e) NULL
  )

  # Write per-outcome instrument files
  sl <- safe_label(label)

  inst_summary_pre  <- pre_harm_summary_bone
  inst_summary_post <- inst_summary %>% dplyr::select(-outcome)

  inst_summary_pre  %>%
    write_csv(paste0("Instrument Metrics - Total Cholesterol for ", sl,
                     " - Pre-Harmonization.csv"))
  inst_summary_post %>%
    write_csv(paste0("Instrument Metrics - Total Cholesterol for ", sl,
                     " - Post-Harmonization.csv"))
  harm_steiger %>%
    mutate(R2.exposure = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
           F.exposure  = (R2.exposure * (samplesize.exposure - 2)) / (1 - R2.exposure)) %>%
    write_csv(paste0("Total Cholesterol Instruments for ", sl, ".csv"))

  list(instruments    = inst_summary,
       mr             = mr_res,
       pleiotropy     = pleio,
       heterogeneity  = het,
       harmonised     = harm_steiger)
}

bone_results <- imap(bone_outcomes, run_bone_mr) %>% purrr::compact()

bone_inst_combined  <- map_dfr(bone_results, "instruments")
bone_mr_combined    <- map_dfr(bone_results, "mr")
bone_pleio_combined <- map_dfr(bone_results, "pleiotropy")
bone_het_combined   <- map_dfr(bone_results, "heterogeneity")
```
:::


### Instrument Strength Across Bone Outcomes


::: {.cell}

```{.r .cell-code}
bone_inst_combined %>%
  kable(caption = paste0(
    "Total cholesterol instruments after harmonisation across bone outcome GWAS. ",
    "Pre-harmonization metrics are identical across outcomes (",
    nrow(instruments.tc.rsid), " rsIDs) and are reported in the per-outcome ",
    "Instrument Metrics CSVs."
  ),
  digits = c(NA, 0, 0, 4, 1, 1, 4, 4, 1))
```

::: {.cell-output-display}


Table: Total cholesterol instruments after harmonisation across bone outcome GWAS. Pre-harmonization metrics are identical across outcomes (360 rsIDs) and are reported in the per-outcome Instrument Metrics CSVs.

|outcome                        | num_snps| samplesize.exposure| cumulative_R2| mean_F| median_F| mean_maf| mean_beta| overall_F|
|:------------------------------|--------:|-------------------:|-------------:|------:|--------:|--------:|---------:|---------:|
|Heel BMD (Morris 2019, UKB)    |      261|              420607|        0.0919|  148.4|     52.7|   0.3499|    0.0305|     163.0|
|Femoral neck BMD (Zheng 2015)  |      262|              420607|        0.0916|  147.4|     52.6|   0.3523|    0.0304|     161.8|
|Fractures (Dönertaş 2021, UKB) |      268|              420607|        0.0939|  147.7|     52.7|   0.3525|    0.0303|     162.6|


:::
:::


### MR Estimates Across Bone Outcomes


::: {.cell}

```{.r .cell-code}
bone_mr_combined %>%
  dplyr::select(-starts_with('id')) %>%
  kable(caption = "MR estimates for total cholesterol on each bone outcome (all methods)",
        digits = c(0, 0, 0, 0, 3, 3, 99))
```

::: {.cell-output-display}


Table: MR estimates for total cholesterol on each bone outcome (all methods)

|outcome                        |exposure                       |method                                                    | nsnp|      b|    se|         pval|
|:------------------------------|:------------------------------|:---------------------------------------------------------|----:|------:|-----:|------------:|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Inverse variance weighted (multiplicative random effects) |  257| -0.051| 0.013| 1.449456e-04|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Inverse variance weighted (fixed effects)                 |  257| -0.051| 0.004| 1.030146e-31|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Robust adjusted profile score (RAPS)                      |  257| -0.041| 0.011| 2.753752e-04|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |MR Egger                                                  |  257| -0.036| 0.021| 8.992305e-02|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Weighted median                                           |  257| -0.026| 0.010| 1.301763e-02|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Weighted mode                                             |  257| -0.027| 0.008| 5.816607e-04|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Inverse variance weighted (multiplicative random effects) |  257| -0.006| 0.020| 7.503772e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Inverse variance weighted (fixed effects)                 |  257| -0.006| 0.018| 7.229353e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Robust adjusted profile score (RAPS)                      |  257| -0.004| 0.021| 8.348923e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |MR Egger                                                  |  257|  0.016| 0.032| 6.221883e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Weighted median                                           |  257| -0.019| 0.031| 5.328331e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Weighted mode                                             |  257|  0.012| 0.029| 6.828918e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Inverse variance weighted (multiplicative random effects) |  264|  0.000| 0.001| 7.405123e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Inverse variance weighted (fixed effects)                 |  264|  0.000| 0.001| 7.293549e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Robust adjusted profile score (RAPS)                      |  264|  0.000| 0.001| 6.952238e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |MR Egger                                                  |  264| -0.001| 0.001| 5.411579e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Weighted median                                           |  264|  0.000| 0.001| 7.700658e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Weighted mode                                             |  264|  0.000| 0.001| 9.076251e-01|


:::
:::


### Primary Estimates Only


::: {.cell}

```{.r .cell-code}
bone_mr_combined %>%
  filter(method == "Inverse variance weighted (multiplicative random effects)") %>%
  dplyr::select(outcome, nsnp, b, se, pval) %>%
  kable(caption = paste0(
    "Primary IVW-RE estimates across bone outcomes — total cholesterol effect ",
    "in standard deviation units (or risk difference per SD for fracture)"
  ),
  digits = c(NA, 0, 4, 4, 6))
```

::: {.cell-output-display}


Table: Primary IVW-RE estimates across bone outcomes — total cholesterol effect in standard deviation units (or risk difference per SD for fracture)

|outcome                        | nsnp|       b|     se|     pval|
|:------------------------------|----:|-------:|------:|--------:|
|Heel BMD (Morris 2019, UKB)    |  257| -0.0509| 0.0134| 0.000145|
|Femoral neck BMD (Zheng 2015)  |  257| -0.0064| 0.0202| 0.750377|
|Fractures (Dönertaş 2021, UKB) |  264|  0.0002| 0.0007| 0.740512|


:::
:::


### Pleiotropy and Heterogeneity Across Bone Outcomes


::: {.cell}

```{.r .cell-code}
bone_pleio_combined %>%
  dplyr::select(-starts_with('id')) %>%
  kable(caption = "MR-Egger intercept tests across bone outcomes")
```

::: {.cell-output-display}


Table: MR-Egger intercept tests across bone outcomes

|outcome                        |exposure                       | egger_intercept|        se|      pval|
|:------------------------------|:------------------------------|---------------:|---------:|---------:|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |      -0.0006185| 0.0007005| 0.3780436|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |      -0.0009350| 0.0010573| 0.3773167|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |       0.0000359| 0.0000342| 0.2950596|


:::

```{.r .cell-code}
bone_het_combined %>%
  dplyr::select(-starts_with('id')) %>%
  kable(caption = "Heterogeneity (Cochran's Q) across bone outcomes",
        digits = c(NA, NA, 0, 3, 0, 4, 1))
```

::: {.cell-output-display}


Table: Heterogeneity (Cochran's Q) across bone outcomes

|outcome                        |exposure                       |method                    |        Q| Q_df| Q_pval|   I2|
|:------------------------------|:------------------------------|:-------------------------|--------:|----:|------:|----:|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |MR Egger                  | 2427.477|  255| 0.0000| 89.5|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Inverse variance weighted | 2434.900|  256| 0.0000| 89.5|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |MR Egger                  |  316.953|  255| 0.0050| 19.5|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Inverse variance weighted |  317.925|  256| 0.0051| 19.5|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |MR Egger                  |  285.846|  262| 0.1490|  8.3|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Inverse variance weighted |  287.047|  263| 0.1476|  8.4|


:::
:::


### Heel BMD — Diagnostic Plots

The Heel BMD (Morris 2019) result is the primary BMD estimate and warrants
detailed diagnostics. Vitamin D was already diagnosed above; the other bone
outcomes (femoral neck BMD, fractures) are secondary or null and are
summarised in the combined results table.


::: {.cell}

```{.r .cell-code}
heelbmd.harmonised <- bone_results[["Heel BMD (Morris 2019, UKB)"]]$harmonised
heelbmd.mr         <- bone_results[["Heel BMD (Morris 2019, UKB)"]]$mr
```
:::


#### Scatter Plot


::: {.cell}

```{.r .cell-code}
ggplot(heelbmd.harmonised, aes(x = beta.exposure, y = beta.outcome)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome,
                    ymax = beta.outcome + 1.96 * se.outcome),
                alpha = 0.5) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96 * se.exposure,
                    xmax = beta.exposure + 1.96 * se.exposure),
                alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 16) +
  labs(x = "Exposure Estimate (Total Cholesterol)",
       y = "Outcome Estimate (Heel BMD)",
       title = "")
```

::: {.cell-output-display}
![](figures/heelbmd-scatter-1.png){width=672}
:::
:::


#### Funnel Plot


::: {.cell}

```{.r .cell-code}
heelbmd.single_snp <- mr_singlesnp(heelbmd.harmonised)

heelbmd.ivw_beta <- heelbmd.mr |>
  filter(method == "Inverse variance weighted (multiplicative random effects)") |>
  pull(b)

heelbmd.y_max <- max(heelbmd.single_snp$se^{-1}, na.rm = TRUE) * 1.1
heelbmd.precision_grid <- seq(0, heelbmd.y_max, length.out = 1000)
heelbmd.bounds_df <- data.frame(
  precision = heelbmd.precision_grid,
  lower     = heelbmd.ivw_beta - 1.96 / heelbmd.precision_grid,
  upper     = heelbmd.ivw_beta + 1.96 / heelbmd.precision_grid
)

ggplot(heelbmd.single_snp, aes(x = b, y = 1/se)) +
  geom_point(size = 1) +
  geom_vline(xintercept = heelbmd.ivw_beta, linetype = "solid",
             color = "#ff7f0e", size = 1) +
  geom_line(data = heelbmd.bounds_df,
            aes(x = lower, y = precision), linetype = "dashed") +
  geom_line(data = heelbmd.bounds_df,
            aes(x = upper, y = precision), linetype = "dashed") +
  labs(x = "Estimate (Beta-IVW)", y = "Precision (1/Standard Error)",
       title = "") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, heelbmd.y_max),
                  xlim = c(min(heelbmd.single_snp$b, na.rm = TRUE),
                           max(heelbmd.single_snp$b, na.rm = TRUE)))
```

::: {.cell-output-display}
![](figures/heelbmd-funnel-1.png){width=672}
:::
:::


#### Leave-One-Out Analysis


::: {.cell}

```{.r .cell-code}
heelbmd.loo_res <- mr_leaveoneout(heelbmd.harmonised)

heelbmd.loo_res |>
  mutate(diff = b - filter(heelbmd.mr,
                           method == "Inverse variance weighted (multiplicative random effects)")$b) |>
  arrange(-abs(diff)) |>
  head() |>
  dplyr::select(SNP, diff, b, se, p) |>
  kable(caption = "Leave-One-Out Results for Heel BMD Analysis (IVW method) for influential SNPs",
        digits = c(0, 5, 5, 5, 5))
```

::: {.cell-output-display}


Table: Leave-One-Out Results for Heel BMD Analysis (IVW method) for influential SNPs

|SNP       |     diff|        b|      se|       p|
|:---------|--------:|--------:|-------:|-------:|
|rs4841132 |  0.00707| -0.04382| 0.01302| 0.00076|
|rs4420638 | -0.00562| -0.05651| 0.01397| 0.00005|
|rs2927472 | -0.00363| -0.05452| 0.01366| 0.00007|
|rs602633  |  0.00349| -0.04740| 0.01370| 0.00054|
|rs174545  |  0.00296| -0.04793| 0.01336| 0.00033|
|rs2737247 |  0.00270| -0.04819| 0.01290| 0.00019|


:::

```{.r .cell-code}
ggplot(heelbmd.loo_res, aes(x = reorder(SNP, -b), y = b)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se),
                width = 0.01, alpha = 0.5) +
  coord_flip() +
  labs(x = "SNP Removed", y = "Estimate (Beta-IVW; leave-one-out)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_text(size = 1))
```

::: {.cell-output-display}
![](figures/heelbmd-loo-1.png){width=672}
:::
:::


Leave-one-out analyses suggested that no individual SNP had a relatively large
influence on the IVW estimate, supporting the robustness of the causal
inference for heel BMD.

#### MR-PRESSO Analysis


::: {.cell}

```{.r .cell-code}
set.seed(2026)
heelbmd.mrpresso <- tryCatch(
  mr_presso(BetaOutcome     = "beta.outcome",
            BetaExposure    = "beta.exposure",
            SdOutcome       = "se.outcome",
            SdExposure      = "se.exposure",
            OUTLIERtest     = TRUE,
            DISTORTIONtest  = TRUE,
            data            = heelbmd.harmonised,
            NbDistribution  = 50000,
            SignifThreshold = 0.05),
  error = function(e) { message("MR-PRESSO failed: ", e$message); NULL }
)

if (!is.null(heelbmd.mrpresso)) {
  heelbmd.mrpresso$`Main MR results` |>
    kable(caption = "MR-PRESSO results for Heel BMD Analysis",
          digits = c(0, 4, 4, 4, 0, 99))

  cat("Global heterogeneity test p-value:",
      heelbmd.mrpresso$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
}
```

::: {.cell-output .cell-output-stdout}

```
Global heterogeneity test p-value: <2e-05 
```


:::

```{.r .cell-code}
heelbmd.mrpresso_rows <- extract_mrpresso_rows(
  heelbmd.mrpresso,
  "Total Cholesterol (UK Biobank)",
  "Heel BMD (Morris 2019, UKB)",
  nrow(heelbmd.harmonised)
)
```
:::


### Forest Plots — All Methods Per Outcome


::: {.cell}

```{.r .cell-code}
ggplot(bone_mr_combined,
       aes(y = method, x = b)) +
  geom_point() +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~outcome, scales = "free_x", ncol = 2) +
  theme_classic(base_size = 12) +
  labs(title = "Total cholesterol effects on bone outcomes",
       y = "", x = "Effect Size (Beta)")
```

::: {.cell-output-display}
![](figures/bone-forest-plots-1.png){width=864}
:::
:::


### Forest Plot — Primary Estimates Across Outcomes


::: {.cell}

```{.r .cell-code}
bone_mr_combined %>%
  filter(method == "Inverse variance weighted (multiplicative random effects)") %>%
  ggplot(aes(y = outcome, x = b)) +
  geom_point(size = 3, colour = color_scheme[1]) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                width = 0.2, colour = color_scheme[1]) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_classic(base_size = 14) +
  labs(title = "Total cholesterol IVW-RE effects on bone outcomes",
       subtitle = "Primary estimate per outcome; 95% CI",
       y = "", x = "Effect Size (Beta)")
```

::: {.cell-output-display}
![](figures/bone-forest-primary-1.png){width=864}
:::
:::


## Summary of Proposed Causal Mechanisms


::: {.cell}

```{.r .cell-code}
# Combine vitamin D with all bone outcomes
all_mr_combined <- bind_rows(vitd.mr, bone_mr_combined) %>%
  filter(method == "Inverse variance weighted (multiplicative random effects)")

ggplot(all_mr_combined, aes(x = outcome, y = b)) +
  geom_point(stat = "identity", size = 3) +
  geom_errorbar(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se), width = 0.2) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic(base_size = 14) +
  labs(title = "Cholesterol effects on calcium-related mechanisms",
       y = "Beta Coefficient", x = "")
```

::: {.cell-output-display}
![](figures/summary-mechanisms-1.png){width=864}
:::
:::


### Primary Mechanisms — Vitamin D and Heel BMD Only

To highlight the two pre-specified mechanistic candidates (vitamin D
synthesis from cholesterol; bone demineralisation), this plot shows only the
two primary outcomes for the calcium hypothesis test.


::: {.cell}

```{.r .cell-code}
primary_mechanisms <- all_mr_combined %>%
  filter(outcome %in% c("Vitamin D (MGI-BioVU LabWAS)",
                        "Heel BMD (Morris 2019, UKB)"))

ggplot(primary_mechanisms, aes(x = outcome, y = b)) +
  geom_point(size = 3, colour = color_scheme[1]) +
  geom_errorbar(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se),
                width = 0.2, colour = color_scheme[1]) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic(base_size = 14) +
  labs(title = "Cholesterol effects on mechanistic outcomes",
       subtitle = "IVW-RE; 95% CI",
       y = "Beta Coefficient (per SD total cholesterol)", x = "")
```

::: {.cell-output-display}
![](figures/summary-primary-mechanisms-1.png){width=768}
:::
:::


## Combined Results — Bone Mechanisms


::: {.cell}

```{.r .cell-code}
# Combine main MR results across all outcomes plus MR-PRESSO rows
all_mr_full <- bind_rows(
  vitd.mr %>% mutate(outcome = "Vitamin D (MGI-BioVU LabWAS)"),
  bone_mr_combined
) %>%
  dplyr::select(-starts_with('id'))

# Append MR-PRESSO rows
all_mrpresso <- bind_rows(
  vitd.mrpresso_rows    %>% mutate(global_pval = as.numeric(global_pval)),
  heelbmd.mrpresso_rows %>% mutate(global_pval = as.numeric(global_pval))
)

# Pleiotropy and heterogeneity rows for the summary file
all_pleio <- bind_rows(
  vitd.pleio %>% mutate(outcome = "Vitamin D (MGI-BioVU LabWAS)"),
  bone_pleio_combined
) %>%
  dplyr::select(-starts_with('id')) %>%
  dplyr::rename(b = egger_intercept) %>%
  mutate(method = "MR-Egger intercept", nsnp = NA_integer_) %>%
  dplyr::select(exposure = exposure, outcome, method, nsnp, b, se, pval)

all_het <- bind_rows(
  vitd.het %>% mutate(outcome = "Vitamin D (MGI-BioVU LabWAS)"),
  bone_het_combined
) %>%
  dplyr::select(-starts_with('id'))

# Master combined results table
mr_results_master <- bind_rows(
  all_mr_full %>%
    mutate(method = fct_recode(as.factor(method),
            "IVW-RE"   = "Inverse variance weighted (multiplicative random effects)",
            "IVW-FE"   = "Inverse variance weighted (fixed effects)",
            "MR-RAPS"  = "Robust adjusted profile score (RAPS)")),
  all_mrpresso %>% dplyr::select(-global_pval)
) %>%
  arrange(outcome, method)

mr_results_master |>
  kable(caption = "MR Results — Bone Mechanisms (all outcomes, all methods)",
        digits = c(0, 0, 0, 0, 4, 4, 99))
```

::: {.cell-output-display}


Table: MR Results — Bone Mechanisms (all outcomes, all methods)

|outcome                        |exposure                       |method                | nsnp|       b|     se|         pval|
|:------------------------------|:------------------------------|:---------------------|----:|-------:|------:|------------:|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |IVW-FE                |  257| -0.0064| 0.0182| 7.229353e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |IVW-RE                |  257| -0.0064| 0.0202| 7.503772e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |MR Egger              |  257|  0.0160| 0.0325| 6.221883e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |MR-RAPS               |  257| -0.0043| 0.0208| 8.348923e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Weighted median       |  257| -0.0192| 0.0308| 5.328331e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Weighted mode         |  257|  0.0118| 0.0290| 6.828918e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |IVW-FE                |  264|  0.0002| 0.0006| 7.293549e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |IVW-RE                |  264|  0.0002| 0.0007| 7.405123e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |MR Egger              |  264| -0.0006| 0.0010| 5.411579e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |MR-RAPS               |  264|  0.0003| 0.0007| 6.952238e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Weighted median       |  264| -0.0003| 0.0011| 7.700658e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Weighted mode         |  264| -0.0001| 0.0010| 9.076251e-01|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |IVW-FE                |  257| -0.0509| 0.0043| 1.030146e-31|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |IVW-RE                |  257| -0.0509| 0.0134| 1.449456e-04|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |MR Egger              |  257| -0.0363| 0.0213| 8.992305e-02|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |MR-PRESSO (Corrected) |  261| -0.0485| 0.0089| 1.277012e-07|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |MR-PRESSO (Raw)       |  261| -0.0520| 0.0133| 1.228511e-04|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |MR-RAPS               |  257| -0.0414| 0.0114| 2.753752e-04|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Weighted median       |  257| -0.0258| 0.0104| 1.301763e-02|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Weighted mode         |  257| -0.0266| 0.0076| 5.816607e-04|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |IVW-FE                |  280| -0.0630| 0.0287| 2.795864e-02|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |IVW-RE                |  280| -0.0630| 0.0326| 5.362791e-02|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |MR Egger              |  280| -0.0700| 0.0533| 1.901778e-01|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |MR-PRESSO (Corrected) |  285|      NA|     NA|           NA|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |MR-PRESSO (Raw)       |  285| -0.0664| 0.0325| 4.225469e-02|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |MR-RAPS               |  280| -0.0672| 0.0336| 4.513899e-02|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |Weighted median       |  280| -0.0051| 0.0512| 9.211422e-01|
|Vitamin D (MGI-BioVU LabWAS)   |Total Cholesterol (UK Biobank) |Weighted mode         |  280| -0.0116| 0.0553| 8.338820e-01|


:::

```{.r .cell-code}
# Write the master results file
mr_results_master %>%
  write_csv("MR Results - Bone Mechanisms.csv")

# Also write the pleiotropy and heterogeneity to companion files
all_pleio %>%
  write_csv("MR Results - Bone Mechanisms - Pleiotropy.csv")
all_het %>%
  write_csv("MR Results - Bone Mechanisms - Heterogeneity.csv")
```
:::


## Hypothesis Testing

Given that we have two hypotheses:

- Cholesterol increases calcium by increasing vitamin D
- Cholesterol increases calcium by decreasing bone mineral density

We performed a Bayesian analysis to determine the posterior probabilities of
the four possible outcomes. For the BMD hypothesis we use the **Morris 2019
heel BMD** result as the primary estimate because it has by far the largest
sample size and therefore the tightest standard error; sensitivity analyses
using the Zheng 2015 femoral neck BMD outcome are reported alongside.


::: {.cell}

```{.r .cell-code}
posterior_prob_direction <- function(beta_hat, se,
                                      direction = c("less", "greater")) {
  direction <- match.arg(direction)
  z <- (0 - beta_hat) / se
  if (direction == "less") {
    return(pnorm(z))
  } else {
    return(1 - pnorm(z))
  }
}

bmd_primary_label <- "Heel BMD (Morris 2019, UKB)"

beta_bmd <- bone_mr_combined %>%
  filter(outcome == bmd_primary_label,
         method == "Inverse variance weighted (multiplicative random effects)") %>%
  pull(b)

se_bmd <- bone_mr_combined %>%
  filter(outcome == bmd_primary_label,
         method == "Inverse variance weighted (multiplicative random effects)") %>%
  pull(se)

beta_vitd <- filter(vitd.mr,
                    method == "Inverse variance weighted (multiplicative random effects)") %>%
  pull(b)
se_vitd <- filter(vitd.mr,
                  method == "Inverse variance weighted (multiplicative random effects)") %>%
  pull(se)

p_h1_true  <- posterior_prob_direction(beta_bmd, se_bmd, "less")
p_h1_false <- 1 - p_h1_true
p_h2_true  <- posterior_prob_direction(beta_vitd, se_vitd, "greater")
p_h2_false <- 1 - p_h2_true

p_both_true    <- p_h1_true * p_h2_true
p_only_h1_true <- p_h1_true * p_h2_false
p_only_h2_true <- p_h1_false * p_h2_true
p_neither_true <- p_h1_false * p_h2_false

outcomes_df <- data.frame(
  Outcome = c("Both True", "Only H1 True", "Only H2 True", "Neither True"),
  Description = c(
    "H1 true (β_BMD < 0) and H2 true (β_VitD > 0)",
    "H1 true (β_BMD < 0) and H2 false (β_VitD <= 0)",
    "H1 false (β_BMD >= 0) and H2 true (β_VitD > 0)",
    "H1 false (β_BMD >= 0) and H2 false (β_VitD <= 0)"
  ),
  Posterior_Probability = c(p_both_true, p_only_h1_true,
                             p_only_h2_true, p_neither_true),
  Percentage = sprintf("%.2f%%",
                       c(p_both_true, p_only_h1_true,
                         p_only_h2_true, p_neither_true) * 100)
)

kable(outcomes_df, format = "simple", digits = 4,
      caption = paste0(
        "Joint posterior probabilities (BMD estimate from ",
        bmd_primary_label, ")"
      ))
```

::: {.cell-output-display}


Table: Joint posterior probabilities (BMD estimate from Heel BMD (Morris 2019, UKB))

Outcome        Description                                         Posterior_Probability  Percentage 
-------------  -------------------------------------------------  ----------------------  -----------
Both True      H1 true (β_BMD < 0) and H2 true (β_VitD > 0)                       0.0268  2.68%      
Only H1 True   H1 true (β_BMD < 0) and H2 false (β_VitD <= 0)                     0.9731  97.31%     
Only H2 True   H1 false (β_BMD >= 0) and H2 true (β_VitD > 0)                     0.0000  0.00%      
Neither True   H1 false (β_BMD >= 0) and H2 false (β_VitD <= 0)                   0.0001  0.01%      


:::
:::


### Bayesian Probability Bar Plot


::: {.cell}

```{.r .cell-code}
bayes_plot_data <- outcomes_df %>%
  mutate(
    Outcome = factor(Outcome,
                     levels = rev(c("Only H1 True", "Both True",
                                    "Only H2 True", "Neither True")))
  )

ggplot(bayes_plot_data,
       aes(x = Posterior_Probability, y = Outcome)) +
  geom_col(fill = color_scheme[1]) +
  geom_text(aes(label = Percentage), hjust = -0.1, size = 4) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, 1)) +
  theme_classic(base_size = 14) +
  labs(title = "Joint posterior probabilities for calcium mechanism hypotheses",
       subtitle = paste0("H1: cholesterol -> reduced BMD; ",
                         "H2: cholesterol -> increased Vitamin D"),
       x = "Posterior Probability", y = "")
```

::: {.cell-output-display}
![](figures/bayes-bar-plot-1.png){width=768}
:::
:::


### Sensitivity — Bayesian Probabilities Using Each BMD Outcome

To assess robustness of the H1 conclusion across different BMD measurement
modalities, we recompute the posterior probabilities using each of the BMD
estimates separately.


::: {.cell}

```{.r .cell-code}
bmd_outcomes_for_bayes <- bone_mr_combined %>%
  filter(method == "Inverse variance weighted (multiplicative random effects)",
         str_detect(outcome, "BMD"))

bayes_sensitivity <- bmd_outcomes_for_bayes %>%
  rowwise() %>%
  mutate(
    p_h1_true     = posterior_prob_direction(b, se, "less"),
    p_only_h1     = p_h1_true * p_h2_false,
    p_both_true   = p_h1_true * p_h2_true,
    `Only H1 (%)` = sprintf("%.2f%%", p_only_h1 * 100),
    `Both H1+H2 (%)` = sprintf("%.2f%%", p_both_true * 100)
  ) %>%
  ungroup() %>%
  dplyr::select(outcome, b, se, pval, p_h1_true,
                `Only H1 (%)`, `Both H1+H2 (%)`)

kable(bayes_sensitivity,
      caption = "Bayesian posterior sensitivity across BMD outcomes",
      digits = c(NA, 4, 4, 4, 4, NA, NA))
```

::: {.cell-output-display}


Table: Bayesian posterior sensitivity across BMD outcomes

|outcome                       |       b|     se|   pval| p_h1_true|Only H1 (%) |Both H1+H2 (%) |
|:-----------------------------|-------:|------:|------:|---------:|:-----------|:--------------|
|Heel BMD (Morris 2019, UKB)   | -0.0509| 0.0134| 0.0001|    0.9999|97.31%      |2.68%          |
|Femoral neck BMD (Zheng 2015) | -0.0064| 0.0202| 0.7504|    0.6248|60.81%      |1.68%          |


:::
:::


### Mathematical Approach

To evaluate the two hypotheses regarding elevated calcium levels:

- H1: lower bone mineral density (BMD, supported by $\beta_{BMD} < 0$)
- H2: higher vitamin D levels (supported by $\beta_{VitD} > 0$)

We applied Bayesian inference using Mendelian randomization (MR) point
estimates and standard errors. The analysis assumed flat (non-informative)
priors on the effect sizes and independence between the effects of BMD and
vitamin D on calcium levels.

For each hypothesis, we modelled the effect size $\beta$ with a flat prior,
$p(\beta) \propto 1$. Given the MR point estimate ($\hat{\beta}$) and standard
error $SE$, the posterior distribution for $\beta$ is approximately Normal:

$$p(\beta | \text{data}) \sim \mathcal{N}(\hat{\beta}, \text{SE}^2)$$

- **H1 (Lower BMD explains elevated calcium)**: The MR estimate from
  Heel BMD (Morris 2019, UKB) is $\hat{\beta}_{\text{BMD}}$ = -0.0509,
  $SE_{BMD}$ = 0.0134. The posterior probability that
  $\beta_{BMD} < 0$ (supporting H1) is 0.9999:

$$P(\beta_{\text{BMD}} < 0 | \text{data}) = \Phi\left(\frac{0 - \hat{\beta}_{\text{BMD}}}{\text{SE}_{\text{BMD}}}\right)$$

- **H2 (Higher Vitamin D)**: The MR estimate is
  $\hat{\beta}_{VitD}$ = -0.063,
  $SE_{VitD}$ = 0.0326. The posterior probability that
  $\beta_{\text{VitD}} > 0$ (supporting H2) is 0.0268:

$$P(\beta_{\text{VitD}} > 0 | \text{data}) = 1 - \Phi\left(\frac{0 - \hat{\beta}_{\text{VitD}}}{\text{SE}_{\text{VitD}}}\right)$$

Thus, the posterior probabilities are estimated at
100% for H1 ($\beta_{\text{BMD}} < 0$) and
2.7% for H2 ($\beta_{VitD} > 0$).

### Joint Posterior Probabilities

Assuming independence between the effects of BMD and vitamin D, we calculated
the joint probabilities for the four possible outcomes:

- **Both True**: 2.68%
- **Only H1 True**: 97.31%
- **Only H2 True**: 0%
- **Neither True**: 0.01%

These probabilities sum to 1, confirming the calculations.

## Summary of MR Analyses

### Robustness of the Heel BMD Effect

The Heel BMD result is robust across all sensitivity analyses. The
MR-Egger intercept was non-significant (intercept =
-6.2\times 10^{-4},
p = 0.38),
indicating no evidence of directional pleiotropy. To address the high
between-SNP heterogeneity (Cochran's Q I² =
89.5%),
we additionally applied MR-PRESSO with NbDistribution = 50,000 to enable
high-precision outlier detection. The global heterogeneity test was
significant (p =
2e-05),
and 45
of 261 SNPs were flagged as pleiotropic outliers.
After removing these outliers, the corrected estimate (β =
-0.0485,
SE = 0.0089,
p = 1.3e-07)
was essentially unchanged from the raw IVW-RE estimate (β =
-0.0509),
with the standard error in fact decreasing because the heterogeneity
contributing to inflated variance was removed. The MR-PRESSO distortion
test, which formally compares raw vs. corrected estimates, returned
p = 0.6,
confirming that outlier-driven pleiotropy is not the source of the
observed cholesterol→BMD effect.

### Vitamin D Does Not Meet the Threshold for a Causal Mediator

In contrast to the heel BMD result, vitamin D failed to meet any of the
pre-specified causal criteria. The IVW-RE point estimate was negative
(β = -0.063,
SE = 0.0326,
p = 0.054),
which is the *opposite* sign from the H2 hypothesis that elevated
cholesterol increases serum calcium via increased vitamin D synthesis.
The robust-method estimates (weighted median
β = -0.0051,
weighted mode β = -0.0116)
attenuated toward zero, suggesting the modestly negative IVW-RE estimate
may itself reflect residual pleiotropy rather than a true effect.
MR-PRESSO confirmed significant global heterogeneity (p =
3e-04),
but at NbDistribution = 50,000 no individual SNPs met the outlier
threshold, indicating diffuse rather than concentrated pleiotropy and
preventing computation of an outlier-corrected estimate. The Bayesian
posterior probability for H2 (cholesterol → increased vitamin D) was
2.7%, which combined with the rejected
direction of effect provides strong evidence against vitamin D as a
mediator of the cholesterol→calcium relationship. We therefore conclude
that vitamin D does not meet the threshold for a causal mediator and
focus the remainder of this analysis on the bone demineralisation
mechanism.

## References

::: {#refs}
:::

## Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.4.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Detroit
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] TwoSampleMR_0.6.29 knitr_1.51         lubridate_1.9.4    forcats_1.0.1     
 [5] stringr_1.6.0      dplyr_1.1.4        purrr_1.2.1        readr_2.1.6       
 [9] tidyr_1.3.2        tibble_3.3.1       ggplot2_4.0.1      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       xfun_0.55          htmlwidgets_1.6.4  psych_2.5.6       
 [5] ggrepel_0.9.6      lattice_0.22-9     tzdb_0.5.0         vctrs_0.6.5       
 [9] tools_4.5.3        generics_0.1.4     curl_7.0.0         parallel_4.5.3    
[13] pkgconfig_2.0.3    Matrix_1.7-4       data.table_1.18.0  RColorBrewer_1.1-3
[17] S7_0.2.1           lifecycle_1.0.5    rootSolve_1.8.2.4  compiler_4.5.3    
[21] farver_2.1.2       mnormt_2.1.1       htmltools_0.5.9    mr.raps_0.4.3     
[25] yaml_2.3.12        pillar_1.11.1      crayon_1.5.3       nlme_3.1-168      
[29] rsnps_0.6.1        tidyselect_1.2.1   digest_0.6.39      nortest_1.0-4     
[33] stringi_1.8.7      labeling_0.4.3     splines_4.5.3      fastmap_1.2.0     
[37] grid_4.5.3         cli_3.6.5          magrittr_2.0.4     dichromat_2.0-0.1 
[41] crul_1.6.0         withr_3.0.2        scales_1.4.0       bit64_4.6.0-1     
[45] timechange_0.3.0   rmarkdown_2.30     bit_4.6.0          otel_0.2.0        
[49] gridExtra_2.3      hms_1.1.4          evaluate_1.0.5     mgcv_1.9-4        
[53] rlang_1.1.7        Rcpp_1.1.1         glue_1.8.0         httpcode_0.3.0    
[57] rstudioapi_0.17.1  vroom_1.6.7        jsonlite_2.0.0     R6_2.6.1          
[61] plyr_1.8.9        
```


:::
:::

