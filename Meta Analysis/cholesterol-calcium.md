---
title: "Meta Analysis of Correlations between Calcium and Cholesterol"
author: "Dave Bridges and Kaelin Loftus"
format:
  html:
    toc: true
    toc-location: left
editor: visual
execute:
  keep-md: true
---





## Data Souces

Located studies from PubMed searches and checking internal references. Manually re-calculated cholesterol to mM when presented in mg/dL


::: {.cell}

```{.r .cell-code}
data.sheet <- 'https://docs.google.com/spreadsheets/d/e/2PACX-1vTyvQnc6bLRLGT6QXEMHxiAQVbK_zag_JIAjvYjTMXINcqdkBwglmg_mlj_k9ml9QsrNQl-tZgy8ACl/pub?gid=1100702568&single=true&output=csv'
library(readr)
data <- read_csv(data.sheet)#from a google sheet
```
:::


The data can be found in the google sheet https://docs.google.com/spreadsheets/d/e/2PACX-1vTyvQnc6bLRLGT6QXEMHxiAQVbK_zag_JIAjvYjTMXINcqdkBwglmg_mlj_k9ml9QsrNQl-tZgy8ACl/pub?gid=1100702568&single=true&output=csv. This script can be found in /Users/davebrid/Documents/GitHub/PrecisionNutrition/Meta Analysis and was most recently run on Thu Apr 20 09:18:10 2023

## Meta Analysis

Analysed data from mean +/- SD of cases and controls


::: {.cell}

```{.r .cell-code}
library(meta)
library(tidyr)
analysis <- metacor(data=data %>% dplyr::filter(!is.na(`r`)),
                   cor = r, 
                   n = n,
                 studlab = Study,
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 hakn = TRUE,
                 title="Calcium and Cholesterol")

forest.meta(analysis,
            fontsize=6,,
            test.overall.random=TRUE)
```

::: {.cell-output-display}
![](figures/ca-chol-meta-1.png){width=672}
:::
:::


# Analysis

There is solid evidence of cross-sectional associations between cholesterol and calcium levels in multiple studies

# Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}
```
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur ... 10.16

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_1.3.0 meta_6.1-0  readr_2.1.3 knitr_1.41 

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0   xfun_0.36          purrr_1.0.1        splines_4.2.2     
 [5] lattice_0.20-45    vctrs_0.5.2        generics_0.1.3     htmltools_0.5.4   
 [9] yaml_2.3.6         utf8_1.2.2         rlang_1.0.6        pillar_1.8.1      
[13] nloptr_2.0.3       DBI_1.1.3          glue_1.6.2         bit64_4.0.5       
[17] lifecycle_1.0.3    stringr_1.5.0      htmlwidgets_1.6.1  evaluate_0.19     
[21] tzdb_0.3.0         fastmap_1.1.0      metafor_3.8-1      parallel_4.2.2    
[25] curl_4.3.3         fansi_1.0.3        Rcpp_1.0.9         vroom_1.6.0       
[29] jsonlite_1.8.4     bit_4.0.5          lme4_1.1-31        hms_1.1.2         
[33] digest_0.6.31      stringi_1.7.12     dplyr_1.0.10       CompQuadForm_1.4.3
[37] grid_4.2.2         mathjaxr_1.6-0     cli_3.6.0          tools_4.2.2       
[41] magrittr_2.0.3     tibble_3.1.8       crayon_1.5.2       pkgconfig_2.0.3   
[45] ellipsis_0.3.2     MASS_7.3-58.1      Matrix_1.5-3       xml2_1.3.3        
[49] assertthat_0.2.1   minqa_1.2.5        rmarkdown_2.19     rstudioapi_0.14   
[53] R6_2.5.1           boot_1.3-28.1      metadat_1.2-0      nlme_3.1-161      
[57] compiler_4.2.2    
```
:::
:::
