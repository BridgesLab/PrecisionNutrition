---
title: "MR Analysis and Validation of SNPs for Calcium"
author: "Dave Bridges"
date: "April 9, 2025"
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
    dev: ["png", "pdf"]  # Remove !expr, just use array syntax
    fig.keep: "all"
execute:
  echo: true
  warning: false
---


::: {.cell}

```{.r .cell-code}
# hide this code chunk
#| echo: false
#| message: false

# defines the se function
se <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

#load these packages, nearly always needed
library(tidyverse)

# sets maize and blue color scheme
color_scheme <- c("#00274c", "#ffcb05")
```
:::


## Purpose

To validate SNPs for calcium GWAS using those identified using UK Biobank.  This script can be found in /Users/davebrid/Documents/GitHub/PrecisionNutrition/Human Genetics and was most recently run on Mon Sep 29 08:44:15 2025

## Data Entry


::: {.cell}

```{.r .cell-code}
instruments.calcium.file <- 'Calcium Instruments from UKBB.csv'
gwas.calcium.file <- 'PheWeb Summary Statistics/phenocode-Ca.tsv.gz'
samplesize.outcome.calcium <- 46100


# loaded and renamed columns
instruments.calcium <- read_csv(instruments.calcium.file) |>
  rename(
    SNP                       = ID,
    id.exposure               = ID,  
    beta.exposure             = BETA,
    se.exposure               = SE,
    effect_allele.exposure    = EA,
    other_allele.exposure     = OA,
    pval.exposure             = P,
    eaf.exposure              = ALT_FREQS,
    samplesize.exposure       = N_exposure
  ) |>
  mutate(exposure="Calcium (UK Biobank)",
         SNP=id.exposure)


gwas.calcium <- read_tsv(gwas.calcium.file) |>
  mutate(ID=paste(chrom, pos, ref,alt, sep=":")) |>
  rename(
    SNP                        = ID,            # or ID if that’s the matching ID
    id.outcome                 = ID,               # can also just set to a string
    beta.outcome               = beta,
    se.outcome                 = sebeta,
    effect_allele.outcome      = alt,   # whichever is effect allele
    other_allele.outcome       = ref,   # whichever is other allele
    pval.outcome               = pval,
    eaf.outcome                = maf,
  ) |>
  mutate(outcome = "Calcium (Michigan GWAS)",
         SNP=id.outcome,# name of your trait)
         samplesize.outcome = samplesize.outcome.calcium)  # sample size for MGI/BioVU for calcium)
```
:::


This presumes the sample sizes was 46100 from Table 1 of https://doi.org/10.1371/journal.pgen.1009077.

Loaded in the instruments for calcium from UK Biobank from the datafile Calcium Instruments from UKBB.csv and the GWAS summary statistics for calcium from the datafile PheWeb Summary Statistics/phenocode-Ca.tsv.gz.


::: {.cell}

```{.r .cell-code}
library(TwoSampleMR)

data <- harmonise_data(instruments.calcium, gwas.calcium, action = 2)
```
:::


Harmonization results

- We used 362 SNPs as instruments for calcium from UK Biobank.
- There were 25 SNPs in common between the exposure and outcome datasets.
- A total of 26 SNPs remained for use after harmonization
- Removed 0 SNPs due to allele mismatches
- Identified 3 palindromic SNPs 

## Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.7

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
 [1] TwoSampleMR_0.6.22 lubridate_1.9.4    forcats_1.0.0      stringr_1.5.2     
 [5] dplyr_1.1.4        purrr_1.1.0        readr_2.1.5        tidyr_1.3.1       
 [9] tibble_3.3.0       ggplot2_4.0.0      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] bit_4.6.0          gtable_0.3.6       jsonlite_2.0.0     crayon_1.5.3      
 [5] compiler_4.5.1     Rcpp_1.1.0         tidyselect_1.2.1   parallel_4.5.1    
 [9] scales_1.4.0       yaml_2.3.10        fastmap_1.2.0      plyr_1.8.9        
[13] R6_2.6.1           generics_0.1.4     knitr_1.50         htmlwidgets_1.6.4 
[17] pillar_1.11.0      RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.6       
[21] stringi_1.8.7      xfun_0.53          S7_0.2.0           bit64_4.6.0-1     
[25] timechange_0.3.0   cli_3.6.5          withr_3.0.2        magrittr_2.0.4    
[29] digest_0.6.37      grid_4.5.1         vroom_1.6.5        rstudioapi_0.17.1 
[33] hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5        data.table_1.17.8 
[37] evaluate_1.0.5     glue_1.8.0         farver_2.1.2       rmarkdown_2.29    
[41] tools_4.5.1        pkgconfig_2.0.3    htmltools_0.5.8.1 
```


:::
:::

