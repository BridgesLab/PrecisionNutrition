---
title: "Analysis of Screening Logs"
author: "Kai-Lin Jen and Dave Bridges"
date: "October 14, 2025"
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
log.file <- 'Screening Log - Screening Log.csv'
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
  Notes = col_character()
))
```
:::


## Analysis

To this point we have screened 

We have decided to review the full text of 177 articles and not to review the full text of 73 articles.  1947 articles have not yet been screened.

### Exclusion Reasons


::: {.cell}

```{.r .cell-code}
screening.log |>
  filter(FullTextView == "No") |>
  group_by(ExclusionReason) |>
  summarize(n = n()) |>
  arrange(desc(n)) |>
  kable()
```

::: {.cell-output-display}


|ExclusionReason                                |  n|
|:----------------------------------------------|--:|
|Not a human study                              | 21|
|Serum cholesterol not mentioned                | 17|
|Case Study                                     | 14|
|Serum calcium not mentioned                    |  8|
|Review article                                 |  6|
|Meta-analysis or systematic review             |  3|
|Retracted                                      |  1|
|Cell Study                                     |  1|
|Neither serum calcium or cholesterol mentioned |  1|
|Systematic review                              |  1|


:::
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

