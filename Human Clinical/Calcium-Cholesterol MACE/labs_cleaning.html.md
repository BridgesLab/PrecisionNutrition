---
title: "Cleaning of Lab Result Data"
author: "Dave Bridges"
date: "2026-03-03"
format:
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
  gfm:
    html-math-method: webtex

theme: journal

execute:
  echo: true
  warning: false

knitr:
  opts_chunk:
    fig-path: "figures/"          # folder for all figure files
    dev: ["png", "pdf"]           # both formats for every plot
    dpi: 300                      # ← this is the correct place for resolution
    dev.args:
      png:
        type: "cairo-png"         # better anti-aliasing / font rendering (optional but recommended)
      pdf:
        family: "sans"            # font family (optional)
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
library(broom)

# sets maize and blue color scheme
color_scheme <- c("#00274c", "#ffcb05")
```
:::




## Purpose

Clean the data for various lab results, for use in other scrupts

## Raw Data




::: {.cell}

```{.r .cell-code}
library(readr) #loads the readr package
labs.filename <- "combined_data/LabResults.csv" #input file(s)

lab.data <- read_csv(labs.filename)
```
:::




These data can be found in /nfs/turbo/precision-health/DataDirect/HUM00268448 - The Interrelationships Between Blood/Cholesterol and Outcomes/2026-03-23 in a file named no file found.  This input file was most recently updated on unknown.  This script was most recently updated on Tue Mar 24 13:46:13 2026.

There are 172531 participants in this dataset with lab values.

## Data Cleaning




::: {.cell}

```{.r .cell-code}
#uncomment to see list of all labs
lab.data |>
  group_by(ORDER_NAME) |>
  count() |>
  arrange(-n) -> lab.summary.table

lab.summary.table |>
  kable(caption="All labs by test name")
```

::: {.cell-output-display}


Table: All labs by test name

|ORDER_NAME                                             |       n|
|:------------------------------------------------------|-------:|
|Estimated GFR, Creatinine-based formula (CKD-EPI 2021) | 1080508|
|PHOSPHORUS LEVEL                                       |  741717|
|25-Hydroxyvitamin D, Total                             |  383188|
|CALCIUM, IONIZED                                       |  107892|
|CREATININE LEVEL                                       |  106230|
|LDL CHOLESTEROL, DIRECT                                |   43302|
|CHOLESTEROL, TOTAL and HDL                             |   28372|
|CREATININE URINE RANDOM                                |   27322|
|CREATININE URINE 24 HOUR                               |   26229|
|CALCIUM URINE 24 HOURS                                 |   21244|
|TRIGLYCERIDES                                          |   19100|
|CALCIUM LEVEL                                          |   14100|
|External Lab EGFR, Creatinine-based formula            |   13289|
|ALBUMIN LEVEL                                          |   10414|
|CHOLESTEROL                                            |    8937|
|Creatinine EX                                          |    8352|
|ALBUMIN, FLUID                                         |    7532|
|POC CREATININE                                         |    6081|
|CALCIUM URINE RANDOM                                   |    3677|
|Calcium EX                                             |    3232|
|TRIGLYCERIDES, FLUID                                   |    2869|
|Procollagen Type I Intact N-Terminal Propeptide (P1NP) |    2132|
|Beta-CrossLaps (B-CTx), Serum                          |    1472|
|Creatinine, 24Hr Urine EX                              |     713|
|NTX-Telopeptide, Urine                                 |     569|
|Triglycerides EX                                       |     540|
|LOW DENSITY LIPOPROTEIN CHOL                           |     539|
|Osteocalcin, Human                                     |     523|
|NTx Serum                                              |     283|
|Ionized Calcium EX                                     |     207|
|HIGH DENSITY LIPOPROTEIN CHOL                          |     205|
|Random Urine Creatinine EX                             |      98|
|Calcium, 24Hr Urine EX                                 |      88|
|Cholesterol EX                                         |      88|
|Voltage-Gated Calcium Channel Ab                       |      40|
|Random Urine Calcium EX                                |      28|
|HDL EX                                                 |      26|
|Phosphate, 24Hr Urine EX                               |      14|
|Parathyroid Hormone, Intact                            |       8|


:::

```{.r .cell-code}
#lab.summary.table |>
  #filter(grepl("ca", ORDER_NAME, ignore.case = TRUE)) #case insensitive has a ca in it
  #filter(grepl("cho", ORDER_NAME, ignore.case = TRUE)) #case insensitive has a ca in it
  #filter(grepl("tri", ORDER_NAME, ignore.case = TRUE)) 

ionized.calcium.labs <- c('CALCIUM, IONIZED','Ionized Calcium EX')
serum.calcium.labs <- c('Calcium EX','CALCIUM LEVEL')
calcium.labs <- c(serum.calcium.labs,ionized.calcium.labs)
tg.labs <- c('TRIGLYCERIDES','TRIGLYCERIDES, FLUID','Triglycerides EX')
ldl.labs <- c('LDL CHOLESTEROL, DIRECT','LOW DENSITY LIPOPROTEIN CHOL')
hdl.labs <- c('HDL EX','HIGH DENSITY LIPOPROTEIN CHOL')
cholesterol.labs <- c('CHOLESTEROL','Cholesterol EX','CHOLESTEROL, TOTAL and HDL')#includes both total and HDL, so duplicate measures, total has to be definitively the higher of the two, they have the same DeID_EncounterID.  need to filter the HDLs into a different category

lab.data.cleaned <-
  lab.data |>
  mutate(VALUE= as.numeric(VALUE)) |> #first make all the values into 
  mutate(test_name = case_when(ORDER_NAME %in% ionized.calcium.labs ~ 'Ionized Calcium',
                               ORDER_NAME %in% serum.calcium.labs ~ 'Serum Calcium',
                               ORDER_NAME %in% tg.labs ~ 'Triglycerides',
                               ORDER_NAME %in% hdl.labs ~ 'HDL-C',
                               ORDER_NAME %in% ldl.labs ~ 'LDL-C',
                               ORDER_NAME %in% cholesterol.labs ~ 'chol_temp' )) %>% #temporary placeholder for cholesterol
  group_by(DeID_EncounterID, test_name) |>
  mutate(test_name = case_when(test_name == 'chol_temp' ~ if_else(VALUE == max(VALUE, na.rm = TRUE),'Total Cholesterol','HDL-C'),
    TRUE ~ test_name
  )) |>
  ungroup() |>
  mutate(value = case_when(UNIT=="mmol/L"&ORDER_NAME %in% calcium.labs ~ VALUE*4.01,
                           UNIT=="MMOL/L"&ORDER_NAME %in% calcium.labs ~ VALUE*4.01,
                          UNIT=="mg/dL"&ORDER_NAME %in% calcium.labs ~ VALUE,
                          ORDER_NAME %in% c(cholesterol.labs,tg.labs,hdl.labs,ldl.labs) ~ VALUE)) |> #always given in mg/dL for calcium
  filter(!is.na(value)) |> #removes tests without values
  filter(!is.na(test_name))  #removes non-cleaned tests

# testing to count how many different units are present now for each test
#lab.data.cleaned |>
#  group_by(test_name, unit_upper = toupper(UNIT)) |>
#  summarise(n = n(), .groups = "drop") |>
#  arrange(test_name, desc(n)) |>
#  kable(caption="Count of the units in each data type")

lab.data.cleaned |>
  group_by(test_name) |>
  summarize(n = length(!is.na(value)),
            mean=mean(value,na.rm=T),
            sd = sd(value,na.rm=T),
            min = min(value,na.rm=T),
            max = max(value,na.rm=T)) |>
  arrange(-n) |>
  kable(caption="Lab results after cleaning and combining")
```

::: {.cell-output-display}


Table: Lab results after cleaning and combining

|test_name         |      n|       mean|          sd|    min|        max|
|:-----------------|------:|----------:|-----------:|------:|----------:|
|Ionized Calcium   | 107810|   4.817992|   0.4850369| 1.3233|    16.1603|
|LDL-C             |  43836| 107.933707|  41.4164068| 0.0000|   524.0000|
|Total Cholesterol |  22107| 188.112589|  50.0432685| 1.0000|   989.0000|
|Triglycerides     |  20959| 225.727134| 451.9433912| 2.0000| 18110.0000|
|HDL-C             |  15465|  62.323375|  33.6794249| 6.0000|   720.0000|
|Serum Calcium     |   6967|   9.549691|   0.8128664| 5.0000|    18.9000|


:::
:::

::: {.cell}

```{.r .cell-code}
encounter.datafile <- 'combined_data/EncounterAll.csv'
encounter.anthro.datafile <- 'combined_data/EncounterAnthropometricsBMI.csv'

#load in encounter metadata
encounter.data <- read_csv(encounter.datafile)
encounter.anthro.data <- read_csv(encounter.anthro.datafile)

merged.lab.data <-
  lab.data.cleaned |>
  left_join(encounter.data, by=c("DeID_PatientID","DeID_EncounterID")) |>
  left_join(encounter.anthro.data, by=c("DeID_PatientID","DeID_EncounterID"))

output.file <- "combined_data/LabResultsCleaned.csv"
merged.lab.data |>
  write_csv(file=output.file)
```
:::




After cleaning are 50091 participants in this dataset with lab values.

The cleaned results are output at **combined_data/LabResultsCleaned.csv**.


## Session Information




::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

Matrix products: default
BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.3/lib64/R/lib/libRblas.so 
LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.3/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Detroit
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] broom_1.0.12    knitr_1.48      lubridate_1.9.3 forcats_1.0.0  
 [5] stringr_1.5.1   dplyr_1.2.0     purrr_1.0.2     readr_2.1.5    
 [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] bit_4.0.5         gtable_0.3.6      jsonlite_1.8.8    crayon_1.5.3     
 [5] compiler_4.4.3    tidyselect_1.2.1  parallel_4.4.3    scales_1.3.0     
 [9] yaml_2.3.9        fastmap_1.2.0     R6_2.5.1          generics_0.1.3   
[13] backports_1.5.0   htmlwidgets_1.6.4 munsell_0.5.1     pillar_1.9.0     
[17] tzdb_0.4.0        rlang_1.1.7       utf8_1.2.4        stringi_1.8.4    
[21] xfun_0.45         bit64_4.0.5       timechange_0.3.0  cli_3.6.3        
[25] withr_3.0.0       magrittr_2.0.3    digest_0.6.36     grid_4.4.3       
[29] vroom_1.6.5       rstudioapi_0.16.0 hms_1.1.3         lifecycle_1.0.5  
[33] vctrs_0.7.1       evaluate_0.24.0   glue_1.8.0        fansi_1.0.6      
[37] colorspace_2.1-0  rmarkdown_2.27    tools_4.4.3       pkgconfig_2.0.3  
[41] htmltools_0.5.8.1
```


:::
:::
