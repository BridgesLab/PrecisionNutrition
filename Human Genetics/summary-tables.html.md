---
title: "Summary Tables and Figures from MR Analyses"
author: "Dave Bridges"
date: "October 17, 2025"
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

To generate summary tables and figures from individual analyses for publicaiton

## Data Entry


::: {.cell}

```{.r .cell-code}
calcium.instrument.summary.file <- 'Instrument Metrics - Calcium - Post-Harmonization (Total Cholesterol).csv' #identical instruments for LDL-Cholesterol outcomes after harmonisation
ldl.instrument.summary.file <- 'Instrument Metrics - LDL Cholesterol - Post-Harmonization.csv'
tc.instrument.summary.file <- 'Instrument Metrics - Total Cholesterol - Post-Harmonization.csv'

bind_rows(read_csv(calcium.instrument.summary.file) |> mutate(Instrument = "Serum Calcium"),
      read_csv(ldl.instrument.summary.file) |> mutate(Instrument = "LDL-Cholesterol"),
      read_csv(tc.instrument.summary.file) |> mutate(Instrument = "Total Cholesterol")) |>
  relocate(Instrument, .before = everything()) |>
  relocate(overall_F, .before= mean_maf) |>
  rename(N=samplesize.exposure,
         SNPs = num_snps,
         `Cumulative R2`=cumulative_R2,
         `Mean F Statistic` = mean_F,
         `Median F Statistic`=median_F,
         `Overall F Statistic`=overall_F,
         `Mean MAF`=mean_maf,
         `Mean Beta Coefficient`=mean_beta) -> instrument.summary

instrument.summary |> 
  kable(caption="Instrument summary after harmonisation")
```

::: {.cell-output-display}


Table: Instrument summary after harmonisation

|Instrument        | SNPs|      N| Cumulative R2| Mean F Statistic| Median F Statistic| Overall F Statistic|  Mean MAF| Mean Beta Coefficient|
|:-----------------|----:|------:|-------------:|----------------:|------------------:|-------------------:|---------:|---------------------:|
|Serum Calcium     |  277| 385066|     0.0638365|          88.8467|           54.59161|            94.72378| 0.3670863|             0.0259874|
|LDL-Cholesterol   |  236| 420607|     0.0886163|         158.3176|           51.14021|           173.19372| 0.3428665|             0.0321334|
|Total Cholesterol |  285| 420607|     0.0972188|         143.7673|           53.49250|           158.81963| 0.3497576|             0.0301266|


:::

```{.r .cell-code}
instrument.summary |>
  mutate(across(ends_with("R2") | ends_with("Coefficient") | ends_with('MAF'), ~ round(., 3))) %>%
  mutate(across(ends_with("Statistic"), ~ round(., 1))) %>%
  # SNPs and N are left unrounded (no action needed) |>
  write_csv("Instrument Metrics - Post-Harmonization.csv")
```
:::


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
 [1] knitr_1.50      lubridate_1.9.4 forcats_1.0.1   stringr_1.5.2  
 [5] dplyr_1.1.4     purrr_1.1.0     readr_2.1.5     tidyr_1.3.1    
 [9] tibble_3.3.0    ggplot2_4.0.0   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] bit_4.6.0          gtable_0.3.6       jsonlite_2.0.0     crayon_1.5.3      
 [5] compiler_4.5.1     tidyselect_1.2.1   parallel_4.5.1     scales_1.4.0      
 [9] yaml_2.3.10        fastmap_1.2.0      R6_2.6.1           generics_0.1.4    
[13] htmlwidgets_1.6.4  pillar_1.11.1      RColorBrewer_1.1-3 tzdb_0.5.0        
[17] rlang_1.1.6        stringi_1.8.7      xfun_0.53          S7_0.2.0          
[21] bit64_4.6.0-1      timechange_0.3.0   cli_3.6.5          withr_3.0.2       
[25] magrittr_2.0.4     digest_0.6.37      grid_4.5.1         vroom_1.6.5       
[29] rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5       
[33] evaluate_1.0.5     glue_1.8.0         farver_2.1.2       rmarkdown_2.29    
[37] tools_4.5.1        pkgconfig_2.0.3    htmltools_0.5.8.1 
```


:::
:::

