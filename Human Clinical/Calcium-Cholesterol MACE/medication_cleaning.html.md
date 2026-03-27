---
title: "Cleaning of Medication Data"
author: "Dave Bridges"
date: "2026-03-25"
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

Clean the data for various medication data, for use in other scrupts

## Raw Data




::: {.cell}

```{.r .cell-code}
library(readr) #loads the readr package
meds.filename <- "combined_data/MedicationOrdersComprehensive.csv" #input file(s)

meds_raw <- read_csv(meds.filename) |>
  mutate(DeID_PatientID = as.character(DeID_PatientID)) 

cat("Total medication rows loaded:", nrow(meds_raw), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Total medication rows loaded: 789424 
```


:::
:::




These data can be found in /nfs/turbo/precision-health/DataDirect/HUM00268448 - The Interrelationships Between Blood/Cholesterol and Outcomes/2026-03-23 in a file named no file found.  This input file was most recently updated on unknown.  This script was most recently updated on Wed Mar 25 14:48:35 2026.

There are 79103 participants in this dataset with meds values.

## Data Cleaning

### Defining Statin Use




::: {.cell}

```{.r .cell-code}
# Match on known statin drug names ONLY — avoids nystatin, bystatin, etc.
statin_pattern <- regex(
  "\\b(atorvastatin|simvastatin|rosuvastatin|pravastatin|
       lovastatin|fluvastatin|pitavastatin|cerivastatin)\\b",
  ignore_case = TRUE
)

meds_statins_raw <- meds_raw |>
  filter(str_detect(MedicationName, statin_pattern))

cat("Statin rows identified:", nrow(meds_statins_raw), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Statin rows identified: 715611 
```


:::

```{.r .cell-code}
cat("Unique patients with any statin order:",
    n_distinct(meds_statins_raw$DeID_PatientID), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Unique patients with any statin order: 68552 
```


:::

```{.r .cell-code}
# Quick check — confirm nystatin not included
cat("Nystatin rows (should be 0):",
    sum(str_detect(meds_statins_raw$MedicationName,
                   regex("nystatin", ignore_case = TRUE))), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Nystatin rows (should be 0): 0 
```


:::
:::




### Parsing Name and Dose




::: {.cell}

```{.r .cell-code}
meds_statins <- meds_statins_raw |>
  mutate(
    # Standardised drug name (lowercase, no brand name parenthetical)
    drug_name = str_extract(MedicationName,
                            statin_pattern) |> str_to_lower(),

    # Dose in mg — extract first number preceding "mg" (case-insensitive)
    dose_mg = str_extract(MedicationName,
                          regex("([0-9]+(?:\\.[0-9]+)?)\\s*mg",
                                ignore_case = TRUE)) |>
      str_extract("[0-9]+(?:\\.[0-9]+)?") |>
      as.numeric(),

    # Parse order date
    order_date = mdy_hm(DeID_OrderStart) |> as_date()
  )

# ============================================================
# BLOCK 1b: Patch dose_mg for the NA rows
# ============================================================

meds_statins <- meds_statins |>
  mutate(
    dose_mg = case_when(

      # Already parsed — keep as-is
      !is.na(dose_mg) ~ dose_mg,

      # Vytorin / ezetimibe-simvastatin combos:
      # "10/40" or "10-40" → second number is the simvastatin dose
      str_detect(MedicationName,
                 regex("ezetimibe.simvastatin|vytorin",
                       ignore_case = TRUE)) ~
        str_extract(MedicationName,
                    regex("(?:[0-9]+[/\\-])([0-9]+)",
                          ignore_case = TRUE)) |>
        str_extract("[0-9]+$") |>
        as.numeric(),

      # Niacin-simvastatin — MinDose is ambiguous (could be niacin dose)
      # Flag as NA and keep a marker column
      str_detect(MedicationName,
                 regex("niacin", ignore_case = TRUE)) ~ NA_real_,

      # Bulk/misc formulations — not classifiable
      str_detect(MedicationName,
                 regex("bulk|misc", ignore_case = TRUE)) ~ NA_real_,

      # Everything else with NA dose: try MinDose
      # (simple "SIMVASTATIN ORAL", salt forms, amlodipine combo, etc.)
      !is.na(MinDose) ~ as.numeric(MinDose),

      # Still NA — leave as NA
      TRUE ~ NA_real_
    ),

    # Flag combos where dose attribution is uncertain
    dose_source = case_when(
      str_detect(MedicationName, regex("ezetimibe.simvastatin|vytorin",
                                       ignore_case = TRUE)) ~ "vytorin_name",
      str_detect(MedicationName, regex("niacin",
                                       ignore_case = TRUE))  ~ "niacin_ambiguous",
      str_detect(MedicationName, regex("amlodipine",
                                       ignore_case = TRUE))  ~ "mindose",
      is.na(str_extract(MedicationName,
                        regex("[0-9]+(?:\\.[0-9]+)?\\s*mg",
                              ignore_case = TRUE)))           ~ "mindose",
      TRUE                                                    ~ "name_parsed"
    )
  )

# Check results
cat("Remaining NA doses after patch:\n")
```

::: {.cell-output .cell-output-stdout}

```
Remaining NA doses after patch:
```


:::

```{.r .cell-code}
meds_statins |>
  filter(is.na(dose_mg)) |>
  count(MedicationName, dose_source, sort = TRUE) |>
  print(n = 20)
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 17 × 3
   MedicationName                           dose_source          n
   <chr>                                    <chr>            <int>
 1 SIMVASTATIN ORAL                         mindose            143
 2 ATORVASTATIN CALCIUM (LIPITOR ORAL)      mindose            141
 3 ATORVASTATIN CALCIUM (ATORVASTATIN ORAL) mindose             93
 4 EZETIMIBE-SIMVASTATIN ORAL               vytorin_name        65
 5 ROSUVASTATIN CALCIUM (CRESTOR ORAL)      mindose             43
 6 PRAVASTATIN SODIUM (PRAVASTATIN ORAL)    mindose             38
 7 SIMVASTATIN (ZOCOR ORAL)                 mindose             23
 8 ROSUVASTATIN CALCIUM (ROSUVASTATIN ORAL) mindose             14
 9 PRAVASTATIN SODIUM (PRAVACHOL ORAL)      mindose              9
10 PITAVASTATIN CALCIUM (LIVALO ORAL)       mindose              8
11 SIMVASTATIN (BULK) MISC                  mindose              5
12 NIACIN-SIMVASTATIN ORAL                  niacin_ambiguous     4
13 simvastatin (bulk) 100 % Powder          mindose              3
14 AMLODIPINE-ATORVASTATIN ORAL             mindose              2
15 AMLODIPINE/ATORVASTATIN (CADUET ORAL)    mindose              1
16 FLUVASTATIN SODIUM (LESCOL ORAL)         mindose              1
17 PITAVASTATIN CALCIUM (PITAVASTATIN ORAL) mindose              1
```


:::

```{.r .cell-code}
cat("\nDose source breakdown:\n")
```

::: {.cell-output .cell-output-stdout}

```

Dose source breakdown:
```


:::

```{.r .cell-code}
meds_statins |> count(dose_source, sort = TRUE) |> print()
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 4 × 2
  dose_source           n
  <chr>             <int>
1 name_parsed      711032
2 mindose            3649
3 vytorin_name        910
4 niacin_ambiguous     20
```


:::

```{.r .cell-code}
cat("\nRe-run intensity classification check (spot-check combos):\n")
```

::: {.cell-output .cell-output-stdout}

```

Re-run intensity classification check (spot-check combos):
```


:::

```{.r .cell-code}
meds_statins |>
  filter(dose_source != "name_parsed") |>
  count(drug_name, dose_mg, dose_source, sort = TRUE) |>
  print(n = 30)
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 59 × 4
   drug_name    dose_mg dose_source          n
   <chr>          <dbl> <chr>            <int>
 1 simvastatin        0 mindose            599
 2 atorvastatin      10 mindose            582
 3 atorvastatin       0 mindose            549
 4 atorvastatin      20 mindose            493
 5 simvastatin       40 vytorin_name       342
 6 atorvastatin      40 mindose            305
 7 simvastatin       20 vytorin_name       242
 8 atorvastatin      NA mindose            237
 9 atorvastatin      80 mindose            198
10 simvastatin       NA mindose            174
11 simvastatin       80 vytorin_name       142
12 simvastatin       10 vytorin_name       119
13 pravastatin        0 mindose            100
14 rosuvastatin       0 mindose             77
15 simvastatin       NA vytorin_name        65
16 rosuvastatin      NA mindose             57
17 pravastatin       NA mindose             47
18 atorvastatin       1 mindose             30
19 simvastatin       20 mindose             28
20 simvastatin       40 mindose             21
21 simvastatin       20 niacin_ambiguous    15
22 pravastatin       40 mindose             12
23 rosuvastatin       5 mindose             12
24 fluvastatin        0 mindose             11
25 rosuvastatin      20 mindose             10
26 simvastatin       10 mindose             10
27 pitavastatin      NA mindose              9
28 rosuvastatin       1 mindose              8
29 rosuvastatin      10 mindose              8
30 pravastatin        1 mindose              7
# ℹ 29 more rows
```


:::
:::




### Quality of Dose Extraction




::: {.cell}

```{.r .cell-code}
cat("\n--- Dose extraction summary ---\n")
```

::: {.cell-output .cell-output-stdout}

```

--- Dose extraction summary ---
```


:::

```{.r .cell-code}
meds_statins |>
  count(drug_name, dose_mg, sort = TRUE) |>
  print(n = 30)
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 59 × 3
   drug_name    dose_mg      n
   <chr>          <dbl>  <int>
 1 atorvastatin      20 148476
 2 atorvastatin      40 148197
 3 atorvastatin      10  75013
 4 atorvastatin      80  68174
 5 simvastatin       20  54312
 6 simvastatin       40  47434
 7 rosuvastatin      10  38260
 8 rosuvastatin       5  24935
 9 rosuvastatin      20  21794
10 pravastatin       20  21218
11 pravastatin       40  18166
12 simvastatin       10  15107
13 rosuvastatin      40  12810
14 pravastatin       10   7440
15 simvastatin       80   4438
16 pravastatin       80   4297
17 simvastatin        5   2397
18 simvastatin        0    599
19 atorvastatin       0    549
20 pitavastatin       1    316
21 fluvastatin       20    243
22 simvastatin       NA    243
23 atorvastatin      NA    237
24 pitavastatin       2    229
25 fluvastatin       40    130
26 pravastatin        0    100
27 pitavastatin       4     88
28 rosuvastatin       0     77
29 fluvastatin       80     71
30 rosuvastatin      NA     57
# ℹ 29 more rows
```


:::

```{.r .cell-code}
cat("\nRows with unparseable dose (NA):",
    sum(is.na(meds_statins$dose_mg)), "\n")
```

::: {.cell-output .cell-output-stdout}

```

Rows with unparseable dose (NA): 594 
```


:::
:::




### Statin Dose Intensity Quantification




::: {.cell}

```{.r .cell-code}
# ACC/AHA 2013 statin intensity categories
# High:     reduces LDL ≥50%
# Moderate: reduces LDL 30-49%
# Low:      reduces LDL <30%

classify_intensity <- function(drug, dose) {
  case_when(
    # High intensity
    drug == "atorvastatin"  & dose >= 40              ~ "high",
    drug == "rosuvastatin"  & dose >= 20              ~ "high",
    # Moderate intensity
    drug == "atorvastatin"  & dose %in% c(10, 20)    ~ "moderate",
    drug == "rosuvastatin"  & dose %in% c(5, 10)     ~ "moderate",
    drug == "simvastatin"   & dose %in% c(20, 40)    ~ "moderate",
    drug == "pravastatin"   & dose >= 40              ~ "moderate",
    drug == "lovastatin"    & dose >= 40              ~ "moderate",
    drug == "fluvastatin"   & dose >= 80              ~ "moderate",
    drug == "pitavastatin"  & dose %in% c(2, 4)      ~ "moderate",
    # Low intensity
    drug == "simvastatin"   & dose <= 10              ~ "low",
    drug == "pravastatin"   & dose <= 20              ~ "low",
    drug == "lovastatin"    & dose <= 20              ~ "low",
    drug == "fluvastatin"   & dose <= 40              ~ "low",
    drug == "pitavastatin"  & dose == 1               ~ "low",
    # Cerivastatin (withdrawn 2001 — rare, treat as moderate if present)
    drug == "cerivastatin"                            ~ "moderate",
    # If dose is NA or doesn't match above, flag for review
    TRUE                                              ~ NA_character_
  )
}

meds_statins <- meds_statins |>
  mutate(
    intensity = classify_intensity(drug_name, dose_mg),
    intensity = factor(intensity, levels = c("low", "moderate", "high"))
  )

cat("Intensity classification:\n")
```

::: {.cell-output .cell-output-stdout}

```
Intensity classification:
```


:::

```{.r .cell-code}
meds_statins |> count(intensity, drug_name, dose_mg) |> print(n = 40)
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 59 × 4
   intensity drug_name    dose_mg      n
   <fct>     <chr>          <dbl>  <int>
 1 low       fluvastatin      0       11
 2 low       fluvastatin     20      243
 3 low       fluvastatin     40      130
 4 low       pitavastatin     1      316
 5 low       pravastatin      0      100
 6 low       pravastatin      1        7
 7 low       pravastatin      5        7
 8 low       pravastatin     10     7440
 9 low       pravastatin     20    21218
10 low       simvastatin      0      599
11 low       simvastatin      1        6
12 low       simvastatin      2        1
13 low       simvastatin      2.5     16
14 low       simvastatin      5     2397
15 low       simvastatin     10    15107
16 moderate  atorvastatin    10    75013
17 moderate  atorvastatin    20   148476
18 moderate  fluvastatin     80       71
19 moderate  pitavastatin     2      229
20 moderate  pitavastatin     4       88
21 moderate  pravastatin     40    18166
22 moderate  pravastatin     80     4297
23 moderate  rosuvastatin     5    24935
24 moderate  rosuvastatin    10    38260
25 moderate  simvastatin     20    54312
26 moderate  simvastatin     40    47434
27 high      atorvastatin    40   148197
28 high      atorvastatin    80    68174
29 high      atorvastatin   100        1
30 high      rosuvastatin    20    21794
31 high      rosuvastatin    30        1
32 high      rosuvastatin    40    12810
33 high      rosuvastatin    80        2
34 high      rosuvastatin   150        1
35 <NA>      atorvastatin     0      549
36 <NA>      atorvastatin     0.5      1
37 <NA>      atorvastatin     1       30
38 <NA>      atorvastatin     2.5      1
39 <NA>      atorvastatin     5       33
40 <NA>      atorvastatin    NA      237
# ℹ 19 more rows
```


:::

```{.r .cell-code}
cat("\nUnclassified (NA intensity):", sum(is.na(meds_statins$intensity)), "\n")
```

::: {.cell-output .cell-output-stdout}

```

Unclassified (NA intensity): 5752 
```


:::

```{.r .cell-code}
# --- Status cleaning ---
# Keep:    Sent, Dispensed, Verified → prescription start events
# End:     Discontinued, Suspend     → prescription stop events
# Exclude: Canceled                  → never filled
# Exclude: Unknown                   → uninterpretable
# Exclude: Completed                 → too rare (n=1304) and ambiguous

meds_statins_clean <- meds_statins |>
  filter(!Status %in% c("Canceled", "Unknown", "Completed")) |>
  mutate(
    event_type = case_when(
      Status %in% c("Sent", "Dispensed", "Verified") ~ "start",
      Status %in% c("Discontinued", "Suspend")       ~ "stop"
    )
  ) |>
  filter(!is.na(order_date)) |>
  arrange(DeID_PatientID, drug_name, order_date)

cat("\nRows after status cleaning:", nrow(meds_statins_clean), "\n")
```

::: {.cell-output .cell-output-stdout}

```

Rows after status cleaning: 623385 
```


:::

```{.r .cell-code}
cat("Start events:", sum(meds_statins_clean$event_type == "start"), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Start events: 523129 
```


:::

```{.r .cell-code}
cat("Stop events: ", sum(meds_statins_clean$event_type == "stop"), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Stop events:  100256 
```


:::
:::




### Patient-Level On/Off Intervals for Statins




::: {.cell}

```{.r .cell-code}
SUPPLY_DAYS <- 90

# Use integer levels for safe comparison inside the loop
INTENSITY_LEVELS <- c("low" = 1L, "moderate" = 2L, "high" = 3L)

build_statin_intervals <- function(df) {

  starts <- df |>
    filter(event_type == "start") |>
    mutate(
      intensity_int = INTENSITY_LEVELS[as.character(intensity)]
    ) |>
    select(order_date, drug_name, intensity_int)

  stops <- df |>
    filter(event_type == "stop") |>
    pull(order_date)

  if (nrow(starts) == 0) return(tibble())

  intervals <- starts |>
    mutate(
      end_provisional = order_date + days(SUPPLY_DAYS),
      stop_date = map(order_date, \(s) {
        valid_stops <- stops[stops > s]
        if (length(valid_stops) == 0) NA_Date_ else min(valid_stops)
      }) |> map_vec(~ if (is.null(.x)) NA_Date_ else .x),
      period_end = pmin(end_provisional,
                        coalesce(stop_date, end_provisional))
    ) |>
    select(period_start = order_date, period_end, intensity_int)

  # Collapse overlapping / abutting intervals
  intervals <- intervals |> arrange(period_start)

  merged    <- vector("list", nrow(intervals))
  current   <- as.list(intervals[1, ])
  merge_idx <- 1L

  for (i in seq_len(nrow(intervals))[-1]) {
    row <- intervals[i, ]
    if (row$period_start <= current$period_end + days(1)) {
      current$period_end    <- max(current$period_end, row$period_end)
      current$intensity_int <- max(current$intensity_int,
                                   row$intensity_int, na.rm = TRUE)
    } else {
      merged[[merge_idx]] <- as_tibble(current)
      merge_idx <- merge_idx + 1L
      current   <- as.list(row)
    }
  }
  merged[[merge_idx]] <- as_tibble(current)

  bind_rows(merged[seq_len(merge_idx)]) |>
    mutate(
      intensity = factor(
        names(INTENSITY_LEVELS)[intensity_int],
        levels = names(INTENSITY_LEVELS)
      )
    ) |>
    select(period_start, period_end, intensity)
}

statin_intervals <- meds_statins_clean |>
  group_by(DeID_PatientID) |>
  group_modify(~ build_statin_intervals(.x)) |>
  ungroup()

cat("Statin interval rows built:", nrow(statin_intervals), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Statin interval rows built: 434960 
```


:::

```{.r .cell-code}
cat("Patients with >=1 statin interval:",
    n_distinct(statin_intervals$DeID_PatientID), "\n")
```

::: {.cell-output .cell-output-stdout}

```
Patients with >=1 statin interval: 65197 
```


:::

```{.r .cell-code}
cat("\nInterval duration summary (days):\n")
```

::: {.cell-output .cell-output-stdout}

```

Interval duration summary (days):
```


:::

```{.r .cell-code}
statin_intervals |>
  mutate(duration = as.numeric(period_end - period_start)) |>
  summarise(
    min    = min(duration),
    median = median(duration),
    mean   = round(mean(duration), 1),
    p90    = quantile(duration, 0.9),
    max    = max(duration)
  ) |> print()
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 1 × 5
    min median  mean   p90   max
  <dbl>  <dbl> <dbl> <dbl> <dbl>
1     1     90  95.4   112  1467
```


:::

```{.r .cell-code}
cat("\nIntensity distribution across intervals:\n")
```

::: {.cell-output .cell-output-stdout}

```

Intensity distribution across intervals:
```


:::

```{.r .cell-code}
statin_intervals |> count(intensity) |> print()
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 4 × 2
  intensity      n
  <fct>      <int>
1 low        26903
2 moderate  255469
3 high      149042
4 <NA>        3546
```


:::
:::

::: {.cell}

```{.r .cell-code}
output.file <- 'combined_data/MedicationOrdersCleanedStatins.csv'
write_csv(statin_intervals, output.file)
```
:::




The cleaned results are output at **combined_data/MedicationOrdersCleanedStatins.csv**.


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
