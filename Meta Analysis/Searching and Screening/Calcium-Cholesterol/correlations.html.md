---
title: "Analysis of Screening Logs"
author: "Kai-Lin Jen, Kaelin Loftus, Cece Occhipinti and Dave Bridges"
date: "May 1, 2026"
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

To quantify screening efforts for the calcium-cholesterol meta-analysis

## Data Entry


::: {.cell}

```{.r .cell-code}
log.file <- 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSTzg4wsY7a4Yz-FNbV74j0_ZTZE9vRmKArv7X7TNwKohKPj1WJRw0EM8aAMuMcfHzrFZbQAZj6nG44/pub?gid=1434844649&single=true&output=csv'
screening.log <- read_csv(log.file,
                          col_types = cols(
  Key = col_character(),
  Title = col_character(),
  DOI = col_character(),
  `Abstract Note` = col_character(),
  DateScreened = col_date(format = "%m/%d/%y"),
  Screener = col_factor(),
  Decision = col_character(),
  FullTextView = col_factor(),
  ExclusionReason = col_factor(),
  Notes = col_character(),
  Correlations = col_factor()
))
```
:::


## Analysis

I wanted to identify papers with the words correlation, association in the title or abstract


::: {.cell}

```{.r .cell-code}
correlation.association.papers <- screening.log %>%
  filter(str_detect(Title, regex("correlation|association", ignore_case = TRUE)) |
           str_detect(`Abstract Note`, regex("correlation|association", ignore_case = TRUE)))

#and papers with just hte word correlation in the title or abstract
correlation.papers <- screening.log %>%
  filter(str_detect(Title, regex("correlation", ignore_case = TRUE)) |
           str_detect(`Abstract Note`, regex("correlation", ignore_case = TRUE)))

correlation.papers |>
  filter(is.na(Correlations)) |>
write_csv("Papers to screen for correlations.csv")
```
:::


There were 649 papers with the words correlation or association in the title or abstract. There were 314 papers with the word correlation in the title or abstract.  Of those 311 remain to be evaluated

## Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Tahoe 26.5.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Detroit
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] knitr_1.51      lubridate_1.9.5 forcats_1.0.1   stringr_1.6.0  
 [5] dplyr_1.2.1     purrr_1.2.2     readr_2.2.0     tidyr_1.3.2    
 [9] tibble_3.3.1    ggplot2_4.0.3   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] bit_4.6.0          gtable_0.3.6       jsonlite_2.0.0     crayon_1.5.3      
 [5] compiler_4.6.0     tidyselect_1.2.1   parallel_4.6.0     scales_1.4.0      
 [9] yaml_2.3.12        fastmap_1.2.0      R6_2.6.1           generics_0.1.4    
[13] curl_7.1.0         htmlwidgets_1.6.4  pillar_1.11.1      RColorBrewer_1.1-3
[17] tzdb_0.5.0         rlang_1.2.0        stringi_1.8.7      xfun_0.57         
[21] S7_0.2.2           bit64_4.8.2        otel_0.2.0         timechange_0.4.0  
[25] cli_3.6.6          withr_3.0.2        magrittr_2.0.5     digest_0.6.39     
[29] grid_4.6.0         vroom_1.7.1        rstudioapi_0.18.0  hms_1.1.4         
[33] lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1        
[37] farver_2.1.2       rmarkdown_2.31     tools_4.6.0        pkgconfig_2.0.3   
[41] htmltools_0.5.9   
```


:::
:::

