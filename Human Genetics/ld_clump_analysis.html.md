---
title: "Calcium and LDLC Clump Analysis from UK Biobank"
author: "Dave Bridges"
date: "2025-06-18"
editor: source
format: 
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
    fig-path: "figures-clumps/"
theme: journal
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


This script analyses the clumped data for the calcium and LDL-C summary statistics from UK Biobank.  The clumps were 

## Data Entry


::: {.cell}

```{.r .cell-code}
calcium.clumps.datafile <- 'biomarkers-30680-both_sexes-irnt-results.clumps'
ldlc.clumps.datafile <- 'biomarkers-30780-both_sexes-irnt-results.clumps'
tc.clumps.datafile <- 'biomarkers-30690-both_sexes-irnt-results.clumps'

calcium.freq.file <- 'biomarkers-30680-all-snps.afreq'
ldlc.freq.file <- 'biomarkers-30780-all-snps.afreq'
tc.freq.file <- 'biomarkers-30690-all-snps.afreq'

calcium.sumstats.file <- 'biomarkers-30680-both_sexes-irnt.tsv'
ldlc.sumstats.file <- 'biomarkers-30780-both_sexes-irnt.tsv'
tc.sumstats.file <- 'biomarkers-30690-both_sexes-irnt.tsv'

calcium.outcome.gwas.file <- 'PheWeb Summary Statistics/phenocode-Ca.tsv.gz'
ldlc.outcome.gwas.file <- 'PheWeb Summary Statistics/phenocode-LDL.tsv.gz'
tc.outcome.gwas.file <- 'PheWeb Summary Statistics/phenocode-Chol.tsv.gz'
```
:::



The calcium summary stats is in biomarkers-30680-both_sexes-irnt.tsv, the total cholesterol in biomarkers-30690-both_sexes-irnt.tsv and the LDL-C summary stats file is in biomarkers-30780-both_sexes-irnt.tsv. 
## LD Clumping

The clumps were generated using the following plink2 options:

```
# or 30780 for LDL-C 
plink2 \
  --pfile 1000G_EUR_merged \
  --clump biomarkers-30680-both_sexes-irnt.tsv \ 
  --clump-snp-field VARIANT \
  --clump-field P \
  --clump-p1 5e-8 \
  --clump-p2 1e-6 \
  --clump-r2 0.01 \
  --clump-kb 10000 \
  --threads 4 \
  --out biomarkers-30680-both_sexes-irnt-results 
```
  
Frequencies we then calculated from these clumped results and are present in the biomarkers-30680-all-snps.afreq, biomarkers-30690-all-snps.afreq and biomarkers-30780-all-snps.afreq
  

## Calcium SNPs


::: {.cell}

```{.r .cell-code}
# Read clumped SNPs and expand out clumped SNPs
calcium.clumps <- read_tsv(calcium.clumps.datafile, col_types = cols()) |>
  rename(P_clumping = P) 

# Read in the outcome GWAS summary statitiscis
calcium.outcome <- read_tsv(calcium.outcome.gwas.file, col_types = cols()) |>
  mutate(ID=paste(chrom, pos, ref,alt, sep=":"))

# First identify the lead SNPs that are in both the clumped SNPs and the outcome GWAS
lead.snps.in.outcome <- intersect(calcium.clumps$ID, calcium.outcome$ID)

# Expanded out clumps from the SP2 column, this is to find SNPs that are buried in the clumps
calcium.clumps.expanded <- calcium.clumps %>%
  separate_rows(SP2, sep = ",")

# Read frequency file (from PLINK --freq output)
calcium.freqs <- read_tsv(calcium.freq.file, col_types = cols())

# Read full summary stats (with BETA)
calcium.sumstats <- read_tsv(calcium.sumstats.file, col_types = cols())

# Found SNPs that matched between clumped SNPs and outcome GWAS
calcium.merged.snps.provisional <-
  left_join(calcium.clumps.expanded, calcium.outcome,
          by = c("SP2"="ID")) |>
  left_join(calcium.freqs, by = c("SP2"="ID","#CHROM"="#CHROM")) |>
  left_join(calcium.sumstats, by = c("SP2"="VARIANT")) |>
  mutate(present_in_outcome = SP2 %in% calcium.outcome$ID) 

# Strategy per clump:
calcium.best_per_clump <- 
  calcium.merged.snps.provisional %>%
  filter(present_in_outcome==T) %>% #lose a lot here
  group_by(ID) %>%
  arrange(
    # 1. prefer if SNP itself is the lead
    desc(SP2 == ID),
    # 2. prefer smaller p-value (stronger in exposure)
    P
  ) %>%
  slice_head(n = 1) %>%   # <-- guarantees one row per group
  ungroup()

# Filter by MAF >= 0.01
calcium.instruments <- 
  calcium.best_per_clump %>%
  filter(ALT_FREQS >= 0.01, ALT_FREQS <= 0.99)  # also exclude rare alleles > 0.99
```
:::


### Instrument Selection Summary for Calcium

| Stage                          | SNPs                          |
|--------------------------------|-------------------------------|
| UKBB SNPs                      | 28987534    |
| After LD Clumping              | 454      |
| After MAF Filtering            | 275 |


::: {.cell}

```{.r .cell-code}
n.calcium <- 385066  

#based on UKBB summary statistics
calcium.instruments <- 
  calcium.instruments %>%
  mutate(
    R2 = 2 * ALT_FREQS * (1 - ALT_FREQS) * BETA^2,
    F = (R2 * (n.calcium - 2)) / (1 - R2)
  )

# Calculate summary metrics
calcium.summary_metrics <- 
  calcium.instruments %>%
  summarise(
    num_snps = n(),
    cumulative_R2 = sum(R2, na.rm = TRUE),
    mean_F = mean(F, na.rm=TRUE),
    median_F = median(F, na.rm=TRUE),
    mean_maf = mean(ALT_FREQS, na.rm = TRUE),
    mean_beta = mean(abs(BETA), na.rm = TRUE),
    overall_F = (cumulative_R2 * (n.calcium - num_snps - 1)) / ((1 - cumulative_R2) * num_snps)
  )

kable(calcium.summary_metrics, caption="Summary of calcium instruments prior to harmonization")
```

::: {.cell-output-display}


Table: Summary of calcium instruments prior to harmonization

| num_snps| cumulative_R2| mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------:|------:|--------:|---------:|---------:|---------:|
|      275|     0.0635624| 89.109| 54.59161| 0.3675692| 0.0260401|  94.97568|


:::

```{.r .cell-code}
calcium.instruments |> 
  rename(EA = alt, ,
         OA = ref,
         CHR = chrom) |>
  mutate(N_exposure=n.calcium) |>
  select(SP2,CHR,POS, EA,OA,BETA,SE,P,ALT_FREQS,N_exposure,R2, `F`) |>
  write_csv("Calcium Instruments from UKBB.csv")
```
:::


For each SNP calculated the $R^2$ and the $F$ statistic using these formulas:

$$
R^2 = 2 \times \text{MAF} \times (1 - \text{MAF}) \times \beta^2
$$

$$
F = \frac{R^2 \times (N - k - 1)}{(1 - R^2) \times k}
$$
Where $MAF$ is the allele frequency from 1000 Genomes, EUR subset, $\beta$ is the effect size of the allele from UK Biobank, $N$ is the sample size of that trait in UK BioBank, and $k$ is the number of instruments (SNPs).  The $F$ statistic is used to assess instrument strength, with values above 10 indicating strong instruments.

## Total Cholesterol SNPs


::: {.cell}

```{.r .cell-code}
# Read clumped SNPs and expand out clumped SNPs
tc.clumps <- read_tsv(tc.clumps.datafile, col_types = cols()) |>
  rename(P_clumping = P) 

# Read in the outcome GWAS summary statitiscis
tc.outcome <- read_tsv(tc.outcome.gwas.file, col_types = cols()) |>
  mutate(ID=paste(chrom, pos, ref,alt, sep=":"))

# First identify the lead SNPs that are in both the clumped SNPs and the outcome GWAS
lead.snps.in.outcome <- intersect(tc.clumps$ID, tc.outcome$ID)

# Expanded out clumps from the SP2 column, this is to find SNPs that are buried in the clumps
tc.clumps.expanded <- tc.clumps %>%
  separate_rows(SP2, sep = ",")

# Read frequency file (from PLINK --freq output)
tc.freqs <- read_tsv(tc.freq.file, col_types = cols())

# Read full summary stats (with BETA)
tc.sumstats <- read_tsv(tc.sumstats.file, col_types = cols())

# Found SNPs that matched between clumped SNPs and outcome GWAS
tc.merged.snps.provisional <-
  left_join(tc.clumps.expanded, tc.outcome,
          by = c("SP2"="ID")) |>
  left_join(tc.freqs, by = c("SP2"="ID","#CHROM"="#CHROM")) |>
  left_join(tc.sumstats, by = c("SP2"="VARIANT")) |>
  mutate(present_in_outcome = SP2 %in% tc.outcome$ID) 

# Strategy per clump:
tc.best_per_clump <- 
  tc.merged.snps.provisional %>%
  filter(present_in_outcome==T) %>% #lose a lot here
  group_by(ID) %>%
  arrange(
    # 1. prefer if SNP itself is the lead
    desc(SP2 == ID),
    # 2. prefer smaller p-value (stronger in exposure)
    P
  ) %>%
  slice_head(n = 1) %>%   # <-- guarantees one row per group
  ungroup()

# Filter by MAF >= 0.01
tc.instruments <- 
  tc.best_per_clump %>%
  filter(ALT_FREQS >= 0.01, ALT_FREQS <= 0.99)  # also exclude rare alleles > 0.99
```
:::


### Instrument Selection Summary for Total Cholesterol

| Stage                          | SNPs                          |
|--------------------------------|-------------------------------|
| UKBB SNPs                      | 28987534    |
| After LD Clumping              | 542      |
| After MAF Filtering            | 277 |


::: {.cell}

```{.r .cell-code}
n.tc <- 420607  

#based on UKBB summary statistics
tc.instruments <- 
  tc.instruments %>%
  mutate(
    R2 = 2 * ALT_FREQS * (1 - ALT_FREQS) * BETA^2,
    F = (R2 * (n.tc - 2)) / (1 - R2)
  )

# Calculate summary metrics
tc.summary_metrics <- 
  tc.instruments %>%
  summarise(
    num_snps = n(),
    cumulative_R2 = sum(R2, na.rm = TRUE),
    mean_F = mean(F, na.rm=TRUE),
    median_F = median(F, na.rm=TRUE),
    mean_maf = mean(ALT_FREQS, na.rm = TRUE),
    mean_beta = mean(abs(BETA), na.rm = TRUE),
    overall_F = (cumulative_R2 * (n.tc - num_snps - 1)) / ((1 - cumulative_R2) * num_snps)
  )

kable(tc.summary_metrics, caption="Summary of total cholesterol instruments prior to harmonization")
```

::: {.cell-output-display}


Table: Summary of total cholesterol instruments prior to harmonization

| num_snps| cumulative_R2|   mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------:|--------:|--------:|---------:|---------:|---------:|
|      277|     0.0954829| 145.2829| 52.68925| 0.3498324| 0.0302356|  160.1837|


:::

```{.r .cell-code}
tc.instruments |> 
  rename(EA = alt, ,
         OA = ref,
         CHR = chrom) |>
  mutate(N_exposure=n.tc) |>
  select(SP2,CHR,POS, EA,OA,BETA,SE,P,ALT_FREQS,N_exposure,R2, `F`) |>
  write_csv("Total Cholesterol Instruments from UKBB.csv")
```
:::


## LDL Cholesterol SNPs


::: {.cell}

```{.r .cell-code}
# Read clumped SNPs and expand out clumped SNPs
ldlc.clumps <- read_tsv(ldlc.clumps.datafile, col_types = cols()) |>
  rename(P_clumping = P) 

# Read in the ouldlcome GWAS summary statitiscis
ldlc.outcome <- read_tsv(ldlc.outcome.gwas.file, col_types = cols()) |>
  mutate(ID=paste(chrom, pos, ref,alt, sep=":"))

# First identify the lead SNPs that are in both the clumped SNPs and the outcome GWAS
lead.snps.in.outcome <- intersect(ldlc.clumps$ID, ldlc.outcome$ID)

# Expanded out clumps from the SP2 column, this is to find SNPs that are buried in the clumps
ldlc.clumps.expanded <- ldlc.clumps %>%
  separate_rows(SP2, sep = ",")

# Read frequency file (from PLINK --freq output)
ldlc.freqs <- read_tsv(ldlc.freq.file, col_types = cols())

# Read full summary stats (with BETA)
ldlc.sumstats <- read_tsv(ldlc.sumstats.file, col_types = cols())

# Found SNPs that matched between clumped SNPs and outcome GWAS
ldlc.merged.snps.provisional <-
  left_join(ldlc.clumps.expanded, ldlc.outcome,
          by = c("SP2"="ID")) |>
  left_join(ldlc.freqs, by = c("SP2"="ID","#CHROM"="#CHROM")) |>
  left_join(ldlc.sumstats, by = c("SP2"="VARIANT")) |>
  mutate(present_in_outcome = SP2 %in% ldlc.outcome$ID) 

# Strategy per clump:
ldlc.best_per_clump <- 
  ldlc.merged.snps.provisional %>%
  filter(present_in_outcome==T) %>% #lose a lot here
  group_by(ID) %>%
  arrange(
    # 1. prefer if SNP itself is the lead
    desc(SP2 == ID),
    # 2. prefer smaller p-value (stronger in exposure)
    P
  ) %>%
  slice_head(n = 1) %>%   # <-- guarantees one row per group
  ungroup()

# Filter by MAF >= 0.01
ldlc.instruments <- 
  ldlc.best_per_clump %>%
  filter(ALT_FREQS >= 0.01, ALT_FREQS <= 0.99)  # also exclude rare alleles > 0.99
```
:::


### Instrument Selection Summary for LDL Cholesterol

| Stage                          | SNPs                          |
|--------------------------------|-------------------------------|
| UKBB SNPs                      | 28987534    |
| After LD Clumping              | 469      |
| After MAF Filtering            | 225 |


::: {.cell}

```{.r .cell-code}
n.ldlc <- 419831  

#based on UKBB summary statistics
ldlc.instruments <- 
  ldlc.instruments %>%
  mutate(
    R2 = 2 * ALT_FREQS * (1 - ALT_FREQS) * BETA^2,
    F = (R2 * (n.tc - 2)) / (1 - R2)
  )

# Calculate summary metrics
ldlc.summary_metrics <- 
  ldlc.instruments %>%
  summarise(
    num_snps = n(),
    cumulative_R2 = sum(R2, na.rm = TRUE),
    mean_F = mean(F, na.rm=TRUE),
    median_F = median(F, na.rm=TRUE),
    mean_maf = mean(ALT_FREQS, na.rm = TRUE),
    mean_beta = mean(abs(BETA), na.rm = TRUE),
    overall_F = (cumulative_R2 * (n.tc - num_snps - 1)) / ((1 - cumulative_R2) * num_snps)
  )

kable(ldlc.summary_metrics, caption="Summary of LDL cholesterol instruments prior to harmonization")
```

::: {.cell-output-display}


Table: Summary of LDL cholesterol instruments prior to harmonization

| num_snps| cumulative_R2|   mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------:|--------:|--------:|---------:|---------:|---------:|
|      225|     0.0813825| 152.4881| 49.95755| 0.3430749|  0.031963|  165.5225|


:::

```{.r .cell-code}
ldlc.instruments |> 
  rename(EA = alt, ,
         OA = ref,
         CHR = chrom) |>
  mutate(N_exposure=n.tc) |>
  select(SP2,CHR,POS, EA,OA,BETA,SE,P,ALT_FREQS,N_exposure,R2, `F`) |>
  write_csv("LDL Cholesterol Instruments from UKBB.csv")
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
