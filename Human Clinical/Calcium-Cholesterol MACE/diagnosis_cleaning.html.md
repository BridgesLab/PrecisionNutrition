---
title: "Cleaning of Diagnoses Data"
author: "Dave Bridges"
date: "2026-03-09"
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

Clean the data for various diagnoses, for use in other scripts

## Raw Data




::: {.cell}

```{.r .cell-code}
library(readr) #loads the readr package
diagnosis.filename <- "combined_data/DiagnosesComprehensiveAll.csv" #input file(s)
diagnosis.data <- read_csv(diagnosis.filename) |>
  mutate(DeID_ProblemObservationDate = mdy_hm(DeID_ProblemObservationDate))

#encounter.datafile <- 'combined_data/EncounterAll.csv'
#encounter.anthro.datafile <- 'combined_data/EncounterAnthropometricsBMI.csv'

#load in encounter metadata
#encounter.data <- read_csv(encounter.datafile)
#encounter.anthro.data <- read_csv(encounter.anthro.datafile)
```
:::




These data can be found in /nfs/turbo/precision-health/DataDirect/HUM00268448 - The Interrelationships Between Blood/Cholesterol and Outcomes/2026-03-23 in a file named no file found.  This input file was most recently updated on unknown.  This script was most recently updated on Tue Mar 24 13:42:25 2026.

There are 70116 participants in this dataset with diagnoses values.

## Data Cleaning

Will first clean the crude diagnosis data and then add in the elixhauser and charson data.

### MACE Data




::: {.cell}

```{.r .cell-code}
icd_pattern_mace <- paste0(
  "^410",          # ICD-9 MI
  "|^433\\.[0-9]1",# ICD-9 Occlusion of precerebral arteries (5th digit = 1)
  "|^434",         # ICD-9 Occlusion of cerebral arteries
  "|^436",         # ICD-9 Acute but ill-defined CVA
  "|^430",         # ICD-9 Subarachnoid haemorrhage
  "|^431",         # ICD-9 Intracerebral haemorrhage
  "|^I21",         # ICD-10 Acute MI
  "|^I22",         # ICD-10 Subsequent MI
  "|^I63",         # ICD-10 Cerebral infarction
  "|^I60",         # ICD-10 Subarachnoid haemorrhage
  "|^I61",         # ICD-10 Intracerebral haemorrhage
  "|^I46"          # ICD-10 Cardiac arrest
)

diagnosis.data.mace <- diagnosis.data %>%
  filter(grepl(icd_pattern_mace, TermCodeMapped, ignore.case = FALSE)) 

mace.diagnoses <-
  diagnosis.data.mace |>
  group_by(DeID_PatientID) |>
  arrange(DeID_ProblemObservationDate) |>
  summarize(MACE.onset = first(DeID_ProblemObservationDate)) |>
  mutate(MACE = TRUE)

library(ggplot2)
ggplot(mace.diagnoses,
       aes(x=MACE.onset)) +
  geom_histogram() +
  theme_classic(base_size=16) +
  labs(x="Year of Incident MACE")
```

::: {.cell-output-display}
![](figures/mace-data-cleaning-1.png){width=2100}
:::
:::





After cleaning are 21661 participants in this dataset with MACE.


### Osteoporosis Data




::: {.cell}

```{.r .cell-code}
icd_pattern_osteoporosis <- paste0(
  "^M80",       # ICD-10 Osteoporosis with pathological fracture
  "|^M81",      # ICD-10 Osteoporosis without pathological fracture
  "|^733\\.0"   # ICD-9 733.0x (covers 733.00–733.09)
)

diagnosis.data.osteoporosis <- diagnosis.data %>%
  filter(grepl(icd_pattern_osteoporosis, TermCodeMapped, ignore.case = FALSE)) 

osteoporosis.diagnoses <-
  diagnosis.data.osteoporosis|>
  group_by(DeID_PatientID) |>
  arrange(DeID_ProblemObservationDate) |>
  summarize(Osteoporosis.onset = first(DeID_ProblemObservationDate)) |>
  mutate(Osteoporosis = TRUE)

library(ggplot2)
ggplot(osteoporosis.diagnoses,
       aes(x=Osteoporosis.onset)) +
  geom_histogram() +
  theme_classic(base_size=16) +
  labs(x="Year of Incident Osteoporosis")
```

::: {.cell-output-display}
![](figures/osteoporosis-data-cleaning-1.png){width=2100}
:::
:::




After cleaning are 22332 participants in this dataset with Osteoporosis

## Writing Out Cleaned Files




::: {.cell}

```{.r .cell-code}
combined.diagnoses <-
  full_join(osteoporosis.diagnoses,mace.diagnoses, by='DeID_PatientID')

cleaned.diagnoses.file <- 'combined_data/DiagnosesCleaned.csv'
write_csv(combined.diagnoses, file=cleaned.diagnoses.file)
```
:::




Wrote this out to combined_data/DiagnosesCleaned.csv

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
 [1] utf8_1.2.4        generics_0.1.3    stringi_1.8.4     hms_1.1.3        
 [5] digest_0.6.36     magrittr_2.0.3    evaluate_0.24.0   grid_4.4.3       
 [9] timechange_0.3.0  fastmap_1.2.0     jsonlite_1.8.8    backports_1.5.0  
[13] fansi_1.0.6       scales_1.3.0      cli_3.6.3         rlang_1.1.7      
[17] crayon_1.5.3      bit64_4.0.5       munsell_0.5.1     withr_3.0.0      
[21] yaml_2.3.9        tools_4.4.3       parallel_4.4.3    tzdb_0.4.0       
[25] colorspace_2.1-0  vctrs_0.7.1       R6_2.5.1          lifecycle_1.0.5  
[29] htmlwidgets_1.6.4 bit_4.0.5         vroom_1.6.5       pkgconfig_2.0.3  
[33] pillar_1.9.0      gtable_0.3.6      glue_1.8.0        xfun_0.45        
[37] tidyselect_1.2.1  rstudioapi_0.16.0 farver_2.1.2      htmltools_0.5.8.1
[41] rmarkdown_2.27    labeling_0.4.3    compiler_4.4.3   
```


:::
:::
