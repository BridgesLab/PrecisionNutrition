---
title: "MR Analyses of Secondary Outcomes for Total Cholesterol on Calcium Homeostasis"
author: "Dave Bridges"
date: "October 16, 2025"
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
library(knitr)

# sets maize and blue color scheme
color_scheme <- c("#00274c", "#ffcb05")
```
:::


## Purpose

To test if SNPs for total cholesterol GWAS identified using UK Biobank relate to other mechanistic or pathological outcomes related to calcium homeostasis and bone health.  This script can be found in /Users/davebrid/Documents/GitHub/PrecisionNutrition/Human Genetics and was most recently run on Sun Oct 19 10:29:13 2025

## Data Entry


::: {.cell}

```{.r .cell-code}
instruments.tc.file <- 'Total Cholesterol Instruments from UKBB.csv'

# loaded and renamed columns
instruments.tc <- read_csv(instruments.tc.file) |>
  rename(
    SNP                       = SP2,
    beta.exposure             = BETA,
    se.exposure               = SE,
    effect_allele.exposure    = EA,
    other_allele.exposure     = OA,
    pval.exposure             = P,
    eaf.exposure              = ALT_FREQS,
    samplesize.exposure       = N_exposure
  ) |>
  mutate(id.exposure="Total Cholesterol (UK Biobank)",
         exposure="Total Cholesterol (UK Biobank)")
```
:::


We used 370 SNPs as instruments for total cholesterol from UK Biobank.  These are found in the Total Cholesterol Instruments from UKBB.csv datafile.

## Mechanistic Outcomes

### Vitamin D Levels

This analysis is to test the hypothesis that total cholesterol impacts 25-hydroxyvitamin D levels, as cholesterol is a precursor for vitamin D synthesis in the skin.  This could positively impact calcium levels indirectly via increased vitamin D.  This was from MGI-BioVU LabWAS data @goldsteinLabWASNovelFindings2020.


::: {.cell}

```{.r .cell-code}
gwas.vitd.file <- 'PheWeb Summary Statistics/phenocode-Vit-D.tsv.gz'
samplesize.outcome.vitd <- 12250  # sample size for MGI/BioVU for calcium from 

gwas.vitd <- read_tsv(gwas.vitd.file) |>
  mutate(ID=paste(chrom, pos, ref,alt, sep=":")) |>
  rename(
    SNP                        = ID,            # or ID if that’s the matching ID
    beta.outcome               = beta,
    se.outcome                 = sebeta,
    effect_allele.outcome      = alt,   # whichever is effect allele
    other_allele.outcome       = ref,   # whichever is other allele
    pval.outcome               = pval,
    eaf.outcome                = maf,
  ) |>
  mutate(id.outcome = "Vitamin D (MGI-BioVU LabWAS)",
         outcome = "Vitamin D (MGI-BioVU LabWAS)",
         samplesize.outcome = samplesize.outcome.vitd)  # sample size for MGI/BioVU for vitamin D)
library(TwoSampleMR)

vitd.data <- harmonise_data(instruments.tc, gwas.vitd, action = 2)
vitd.data_steiger <- steiger_filtering(vitd.data)
vitd.mr <- mr(vitd.data_steiger,
                         method_list = c("mr_ivw", 
                                         "mr_egger_regression",
                                         "mr_weighted_median", 
                                         "mr_weighted_mode"))

vitd.mr |> select(-starts_with('id')) |> 
  kable(caption="MR Results for Total Cholesterol - Vitamin D Analysis",
        digits=c(0,0,0,0,3,3,99))
```

::: {.cell-output-display}


Table: MR Results for Total Cholesterol - Vitamin D Analysis

|outcome                      |exposure                       |method                    | nsnp|      b|    se|       pval|
|:----------------------------|:------------------------------|:-------------------------|----:|------:|-----:|----------:|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Inverse variance weighted |  280| -0.063| 0.033| 0.05362791|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |MR Egger                  |  280| -0.070| 0.053| 0.19017784|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted median           |  280| -0.005| 0.049| 0.91795939|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted mode             |  280| -0.012| 0.054| 0.83022348|


:::

```{.r .cell-code}
ggplot(vitd.mr, aes(y=method,x=b)) +
  geom_point() +
  geom_errorbar(aes(xmin=b-1.96*se, xmax=b+1.96*se), width=0.2) +
  theme_classic(base_size=16) +
  labs(title="25-Hydroxyvitamin D Levels (LabWAS)",
       y="",
       x="Effect Size (Beta)") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") 
```

::: {.cell-output-display}
![](figures/mr-vitd-1.png){width=672}
:::
:::


## Pathological Outcomes

### Bone Mineral Density

#### GEFOS Consortium Lumbar Spine BMD

Lumbar spine tends to reflect trabecular bone, which is metabolically active, hormonally responsive, and more sensitive to lipid and endocrine perturbations. Multiple epidemiologic studies report inverse associations between total cholesterol and LS-BMD, particularly among older adults and postmenopausal women @huAssociationsSerumTotal2023 and @fangNegativeAssociationTotal2022.  We used the pooled (not sex-specific) estimates from the GEFOS consortium's pooled meta-analysis @estradaGenomewideMetaanalysisIdentifies2012.

### Fracture Risk

## Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.7.1

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
 [1] TwoSampleMR_0.6.22 knitr_1.50         lubridate_1.9.4    forcats_1.0.1     
 [5] stringr_1.5.2      dplyr_1.1.4        purrr_1.1.0        readr_2.1.5       
 [9] tidyr_1.3.1        tibble_3.3.0       ggplot2_4.0.0      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] generics_0.1.4     lattice_0.22-7     stringi_1.8.7      hms_1.1.3         
 [5] digest_0.6.37      magrittr_2.0.4     evaluate_1.0.5     grid_4.5.1        
 [9] timechange_0.3.0   RColorBrewer_1.1-3 fastmap_1.2.0      plyr_1.8.9        
[13] jsonlite_2.0.0     scales_1.4.0       mnormt_2.1.1       cli_3.6.5         
[17] rlang_1.1.6        crayon_1.5.3       bit64_4.6.0-1      withr_3.0.2       
[21] yaml_2.3.10        tools_4.5.1        parallel_4.5.1     tzdb_0.5.0        
[25] vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.4    htmlwidgets_1.6.4 
[29] bit_4.6.0          vroom_1.6.5        psych_2.5.6        pkgconfig_2.0.3   
[33] pillar_1.11.1      gtable_0.3.6       glue_1.8.0         data.table_1.17.8 
[37] Rcpp_1.1.0         xfun_0.53          tidyselect_1.2.1   rstudioapi_0.17.1 
[41] farver_2.1.2       nlme_3.1-168       htmltools_0.5.8.1  labeling_0.4.3    
[45] rmarkdown_2.29     compiler_4.5.1     S7_0.2.0          
```


:::
:::

