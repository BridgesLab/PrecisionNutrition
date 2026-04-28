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
run on Mon Apr 27 20:31:46 2026.

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
all_rsids <- bind_rows(rsid_list)

# Ensure types match before join
all_rsids <- all_rsids %>%
  mutate(CHR = as.integer(CHR), BP = as.integer(BP))

tc.positions <- tc.positions %>%
  mutate(CHR = as.integer(CHR), BP = as.integer(BP))

# Now do the join with explicit type safety
instruments.tc.rsid <- instruments.tc %>%
  separate(SNP, into = c("CHR", "BP", "ALLELE0", "ALLELE1"),
           sep = ":", remove = TRUE) %>%
  mutate(CHR = as.integer(CHR), BP = as.integer(BP)) %>%
  left_join(all_rsids, by = c("CHR", "BP")) %>%
  filter(!is.na(RSID)) %>%
  select(
    effect_allele.exposure,
    other_allele.exposure,
    beta.exposure,
    se.exposure,
    pval.exposure,
    eaf.exposure,
    samplesize.exposure,
    id.exposure,
    exposure,
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

```{.r .cell-code}
cat("Loss:                              ",
    nrow(instruments.tc) - nrow(instruments.tc.rsid),
    "(", round(100 * (1 - nrow(instruments.tc.rsid) / nrow(instruments.tc)), 1),
    "%)\n", sep = "")
```

::: {.cell-output .cell-output-stdout}

```
Loss:                              10(2.7%)
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

# Instrument strength
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

kable(vitd.exposure.summary,
      caption = "Summary of total cholesterol instruments after harmonisation for vitamin D analysis")
```

::: {.cell-output-display}


Table: Summary of total cholesterol instruments after harmonisation for vitamin D analysis

| num_snps| samplesize.exposure| cumulative_R2|   mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------------:|-------------:|--------:|--------:|---------:|---------:|---------:|
|      285|              420607|     0.0972188| 143.7673|  53.4925| 0.3497576| 0.0301266|  158.8196|


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
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted median                                           |  280| -0.005| 0.051| 0.92162624|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted mode                                             |  280| -0.012| 0.051| 0.82031569|


:::

```{.r .cell-code}
mr_pleiotropy_test(vitd.data_steiger) |>
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
mr_heterogeneity(vitd.data_steiger) |>
  dplyr::select(-starts_with('id')) |>
  mutate(I2 = pmax(0, (Q - Q_df) / Q) * 100) |>
  kable(caption = "MR Heterogeneity Results for Total Cholesterol - Vitamin D Analysis",
        digits = c(0, 0, 0, 3, 3, 99))
```

::: {.cell-output-display}


Table: MR Heterogeneity Results for Total Cholesterol - Vitamin D Analysis

|outcome                      |exposure                       |method                    |       Q| Q_df|       Q_pval| I2|
|:----------------------------|:------------------------------|:-------------------------|-------:|----:|------------:|--:|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |MR Egger                  | 361.847|  278| 0.0005223738| 23|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Inverse variance weighted | 361.883|  279| 0.0005999283| 23|


:::

```{.r .cell-code}
ggplot(vitd.mr, aes(y = method, x = b)) +
  geom_point() +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0.2) +
  theme_classic(base_size = 16) +
  labs(title = "25-OH Vitamin D (LabWAS)", y = "", x = "Effect Size (Beta)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")
```

::: {.cell-output-display}
![](figures/mr-vitd-1.png){width=672}
:::
:::


## Pathological Outcomes

### Bone Outcome GWAS Selection

The previous analysis used local PheWeb summary statistics from older GEFOS
releases (Estrada 2012; Zheng 2015). These files had three problems: small
sample sizes (32k–53k participants), build-37 marker IDs requiring manual
rsID conversion, and dramatic instrument loss during harmonisation (only 23%
of instruments retained for the 2015 LS-BMD analysis). To address these
issues, this version queries three more recent BMD and fracture GWAS
directly via the OpenGWAS API, using LD proxies (r²>0.8) to recover SNPs
that are not directly present in each outcome dataset.

The three bone-related outcomes are treated as primary analyses because
they ask different but complementary questions:

- **Heel BMD (Morris 2019, UKB, n=426,824)**: highest-powered BMD GWAS
  available, derived from quantitative ultrasound; cortical-bone-dominated
- **Femoral neck and lumbar spine BMD (Estrada 2012, GEFOS)**: DXA-based,
  trabecular-bone-dominated, lower powered but more biologically specific
- **Fracture risk (Trajanoska 2018, GEFOS)**: clinical endpoint integrating
  bone density, geometry, and fall propensity

We retain Estrada 2012 alongside Morris 2019 because heel ultrasound BMD
and DXA-BMD measure different aspects of bone biology — heel BMD reflects
predominantly cortical bone (density and microarchitecture of the calcaneus),
whereas lumbar and femoral DXA-BMD reflect a mix of cortical and trabecular
compartments more relevant to age-related demineralisation. If both point in
the same direction, the convergence is strong evidence; if they diverge,
the result is informative about which bone compartment is affected.

### Querying OpenGWAS for BMD and Fracture Outcomes


::: {.cell}

```{.r .cell-code}
library(ieugwasr)

# OpenGWAS IDs for BMD and fracture outcomes
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
# Fetch outcome data via OpenGWAS for each GWAS
fetch_bone_outcome <- function(gwas_id, label) {
  cat("  Fetching", label, "(", gwas_id, ")...\n")
  
  outcome_dat <- tryCatch(
    extract_outcome_data(
      snps     = instruments.tc.rsid$SNP,
      outcomes = gwas_id,
      proxies  = TRUE,
      rsq      = 0.8,
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
# Coverage summary
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
# Use rsID-keyed instruments for bone analyses (OpenGWAS returns rsIDs)
instruments.tc.for.bone <- instruments.tc.rsid

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
  
  list(instruments = inst_summary,
       mr          = mr_res,
       pleiotropy  = pleio,
       heterogeneity = het,
       harmonised  = harm_steiger)
}

bone_results <- imap(bone_outcomes, run_bone_mr) %>% purrr::compact()

# Combine across outcomes
bone_inst_combined <- map_dfr(bone_results, "instruments")
bone_mr_combined   <- map_dfr(bone_results, "mr")
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
    "Compared to the prior local-file analysis (84 SNPs for LS-BMD-2015 / ",
    "181 for fractures), the OpenGWAS approach with LD proxies recovers ",
    "substantially more instruments."
  ),
  digits = c(NA, 0, 0, 4, 1, 1, 4, 4, 1))
```

::: {.cell-output-display}


Table: Total cholesterol instruments after harmonisation across bone outcome GWAS. Compared to the prior local-file analysis (84 SNPs for LS-BMD-2015 / 181 for fractures), the OpenGWAS approach with LD proxies recovers substantially more instruments.

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
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Weighted median                                           |  257| -0.026| 0.010| 9.624763e-03|
|Heel BMD (Morris 2019, UKB)    |Total Cholesterol (UK Biobank) |Weighted mode                                             |  257| -0.027| 0.007| 4.086711e-04|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Inverse variance weighted (multiplicative random effects) |  257| -0.006| 0.020| 7.503772e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Inverse variance weighted (fixed effects)                 |  257| -0.006| 0.018| 7.229353e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Robust adjusted profile score (RAPS)                      |  257| -0.004| 0.021| 8.348923e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |MR Egger                                                  |  257|  0.016| 0.032| 6.221883e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Weighted median                                           |  257| -0.019| 0.031| 5.420855e-01|
|Femoral neck BMD (Zheng 2015)  |Total Cholesterol (UK Biobank) |Weighted mode                                             |  257|  0.012| 0.030| 6.915255e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Inverse variance weighted (multiplicative random effects) |  264|  0.000| 0.001| 7.405123e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Inverse variance weighted (fixed effects)                 |  264|  0.000| 0.001| 7.293549e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Robust adjusted profile score (RAPS)                      |  264|  0.000| 0.001| 6.951942e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |MR Egger                                                  |  264| -0.001| 0.001| 5.411579e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Weighted median                                           |  264|  0.000| 0.001| 7.659518e-01|
|Fractures (Dönertaş 2021, UKB) |Total Cholesterol (UK Biobank) |Weighted mode                                             |  264|  0.000| 0.001| 9.123676e-01|


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
    "in standard deviation units (or log odds for fracture)"
  ),
  digits = c(NA, 0, 4, 4, 6))
```

::: {.cell-output-display}


Table: Primary IVW-RE estimates across bone outcomes — total cholesterol effect in standard deviation units (or log odds for fracture)

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


## Hypothesis Testing

Given that we have two hypotheses:

- Cholesterol increases calcium by increasing vitamin D
- Cholesterol increases calcium by decreasing bone mineral density

We performed a Bayesian analysis to determine the posterior probabilities of
the four possible outcomes. For the BMD hypothesis we use the **Morris 2019
heel BMD** result as the primary estimate because it has by far the largest
sample size and therefore the tightest standard error; sensitivity analyses
using the Estrada 2012 DXA-BMD outcomes are reported alongside.


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

# Use Morris 2019 heel BMD as primary BMD estimate (largest sample size)
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

# Individual posterior probabilities
p_h1_true  <- posterior_prob_direction(beta_bmd, se_bmd, "less")
p_h1_false <- 1 - p_h1_true
p_h2_true  <- posterior_prob_direction(beta_vitd, se_vitd, "greater")
p_h2_false <- 1 - p_h2_true

# Joint probabilities
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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ieugwasr_1.1.0                          
 [2] SNPlocs.Hsapiens.dbSNP144.GRCh37_0.99.20
 [3] BSgenome_1.76.0                         
 [4] rtracklayer_1.68.0                      
 [5] BiocIO_1.18.0                           
 [6] Biostrings_2.76.0                       
 [7] XVector_0.48.0                          
 [8] GenomicRanges_1.60.0                    
 [9] GenomeInfoDb_1.44.3                     
[10] IRanges_2.42.0                          
[11] S4Vectors_0.46.0                        
[12] BiocGenerics_0.54.1                     
[13] generics_0.1.4                          
[14] TwoSampleMR_0.6.29                      
[15] knitr_1.51                              
[16] lubridate_1.9.4                         
[17] forcats_1.0.1                           
[18] stringr_1.6.0                           
[19] dplyr_1.1.4                             
[20] purrr_1.2.1                             
[21] readr_2.1.6                             
[22] tidyr_1.3.2                             
[23] tibble_3.3.1                            
[24] ggplot2_4.0.1                           
[25] tidyverse_2.0.0                         

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1            psych_2.5.6                
 [3] rootSolve_1.8.2.4           farver_2.1.2               
 [5] S7_0.2.1                    bitops_1.0-9               
 [7] fastmap_1.2.0               RCurl_1.98-1.17            
 [9] GenomicAlignments_1.44.0    XML_3.99-0.20              
[11] digest_0.6.39               timechange_0.3.0           
[13] lifecycle_1.0.5             magrittr_2.0.4             
[15] compiler_4.5.3              rlang_1.1.7                
[17] tools_4.5.3                 yaml_2.3.12                
[19] data.table_1.18.0           S4Arrays_1.8.1             
[21] labeling_0.4.3              mr.raps_0.4.3              
[23] htmlwidgets_1.6.4           bit_4.6.0                  
[25] mnormt_2.1.1                curl_7.0.0                 
[27] DelayedArray_0.34.1         plyr_1.8.9                 
[29] RColorBrewer_1.1-3          abind_1.4-8                
[31] BiocParallel_1.42.2         httpcode_0.3.0             
[33] withr_3.0.2                 grid_4.5.3                 
[35] scales_1.4.0                SummarizedExperiment_1.38.1
[37] crul_1.6.0                  dichromat_2.0-0.1          
[39] cli_3.6.5                   rmarkdown_2.30             
[41] crayon_1.5.3                otel_0.2.0                 
[43] rstudioapi_0.17.1           httr_1.4.7                 
[45] tzdb_0.5.0                  rjson_0.2.23               
[47] rsnps_0.6.1                 splines_4.5.3              
[49] parallel_4.5.3              restfulr_0.0.16            
[51] matrixStats_1.5.0           vctrs_0.6.5                
[53] Matrix_1.7-4                jsonlite_2.0.0             
[55] hms_1.1.4                   bit64_4.6.0-1              
[57] ggrepel_0.9.6               nortest_1.0-4              
[59] glue_1.8.0                  codetools_0.2-20           
[61] stringi_1.8.7               gtable_0.3.6               
[63] UCSC.utils_1.4.0            pillar_1.11.1              
[65] htmltools_0.5.9             GenomeInfoDbData_1.2.14    
[67] R6_2.6.1                    Biobase_2.68.0             
[69] vroom_1.6.7                 evaluate_1.0.5             
[71] lattice_0.22-9              Rsamtools_2.24.1           
[73] Rcpp_1.1.1                  SparseArray_1.8.1          
[75] gridExtra_2.3               nlme_3.1-168               
[77] xfun_0.55                   MatrixGenerics_1.20.0      
[79] pkgconfig_2.0.3            
```


:::
:::

