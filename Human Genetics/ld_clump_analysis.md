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

calcium.freq.file <- 'biomarkers-30680-both_sexes-irnt-results-clumped-snps.afreq'
ldlc.freq.file <- 'biomarkers-30780-both_sexes-irnt-results-clumped-snps.afreq'

calcium.sumstats.file <- 'biomarkers-30680-both_sexes-irnt.tsv'
ldlc.sumstats.file <- 'biomarkers-30780-both_sexes-irnt.tsv'
```
:::


The calcium summary stats is in biomarkers-30680-both_sexes-irnt.tsv and the LDL-C summary stats file is in biomarkers-30780-both_sexes-irnt.tsv.  The clumps were generated using the following plink2 options:

  --clump-id-field VARIANT
  --clump-kb 10000
  --clump-p-field P
  --clump-p1 5e-8
  --clump-p2 1e-6
  --clump-r2 0.001
  --out biomarkers-30780-both_sexes-irnt-results
  --pfile 1000G_EUR_merged
  --threads 4
  
Frequencies we then calculated from these clumped results and are present in the biomarkers-30680-both_sexes-irnt-results-clumped-snps.afreq and biomarkers-30780-both_sexes-irnt-results-clumped-snps.afreq
  

## Calcium SNPs


::: {.cell}

```{.r .cell-code}
# Read clumped SNPs
calcium.clumps <- read_tsv(calcium.clumps.datafile, col_types = cols())

# Read frequency file (from PLINK --freq output)
calcium.freqs <- read_tsv(calcium.freq.file, col_types = cols())

# Read full summary stats (with BETA)
calcium.sumstats <- read_tsv(calcium.sumstats.file, col_types = cols())

# Inspect the column names to confirm merging keys and data availability
# head(calcium.clumps)
# head(calcium.freqs)
# head(calcium.sumstats)

# The clump file 'SNP' column should match frequency 'ID' column
# The sumstats 'VARIANT' column should match clumps 'SNP' column

# Merge clumps with frequency (to get MAF)
calcium.clumps_freq <- calcium.clumps %>%
  inner_join(calcium.freqs %>% select(ID, ALT_FREQS), by = "ID")

# Merge with summary stats (to get beta)
calcium.clumps_freq_beta <- calcium.clumps_freq %>%
  inner_join(calcium.sumstats %>% rename(ID = VARIANT) %>% select(ID, BETA), by = "ID")

# Filter by MAF >= 0.01
calcium.instruments <- calcium.clumps_freq_beta %>%
  filter(ALT_FREQS >= 0.01, ALT_FREQS <= 0.99)  # also exclude rare alleles > 0.99

n.calcium <- 431888  # replace with your actual exposure GWAS sample size

calcium.instruments <- calcium.instruments %>%
  mutate(
    R2 = 2 * ALT_FREQS * (1 - ALT_FREQS) * BETA^2,
    F = (R2 * (n.calcium - 2)) / (1 - R2)
  )

# Calculate summary metrics
calcium.summary_metrics <- calcium.instruments %>%
  summarise(
    num_snps = n(),
    cumulative_R2 = sum(R2, na.rm = TRUE),
    mean_F = mean(F, na.rm=TRUE),
    median_F = median(F, na.rm=TRUE),
    mean_maf = mean(ALT_FREQS, na.rm = TRUE),
    mean_beta = mean(abs(BETA), na.rm = TRUE),
    overall_F = (cumulative_R2 * (n.calcium - num_snps - 1)) / ((1 - cumulative_R2) * num_snps)
  )

kable(calcium.summary_metrics, caption="Summary of calcium instruments")
```

::: {.cell-output-display}
Table: Summary of calcium instruments

| num_snps| cumulative_R2|   mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------:|--------:|--------:|---------:|---------:|---------:|
|      362|     0.0963254| 115.0375| 69.57594| 0.3208676| 0.0328442|   127.065|
:::
:::


For each SNP calculated the $R^2$ and the $F$ statistic using these formulas:

$$
R^2 = 2 \times \text{MAF} \times (1 - \text{MAF}) \times \beta^2
$$

$$
F = \frac{R^2 \times (N - k - 1)}{(1 - R^2) \times k}
$$
Where $MAF$ is the allele frequency from 1000 Genomes, EUR subset, $\beta$ is the effect size of the allele from UK Biobank, $N$ is the sample size of that trait in UK BioBank, and $k$ is the number of instruments (SNPs).  The $F$ statistic is used to assess instrument strength, with values above 10 indicating strong instruments.

## LDL-C SNPs


::: {.cell}

```{.r .cell-code}
# Read clumped SNPs
ldlc.clumps <- read_tsv(ldlc.clumps.datafile, col_types = cols())

# Read frequency file (from PLINK --freq output)
ldlc.freqs <- read_tsv(ldlc.freq.file, col_types = cols())

# Read full summary stats (with BETA)
ldlc.sumstats <- read_tsv(ldlc.sumstats.file, col_types = cols())

# Inspect the column names to confirm merging keys and data availability
# head(ldlc.clumps)
# head(ldlc.freq)
# head(ldlc.sumstats)

# The clump file 'SNP' column should match frequency 'ID' column
# The sumstats 'VARIANT' column should match clumps 'SNP' column

# Merge clumps with frequency (to get MAF)
ldlc.clumps_freq <- ldlc.clumps %>%
  inner_join(ldlc.freqs %>% select(ID, ALT_FREQS), by = "ID")

# Merge with summary stats (to get beta)
ldlc.clumps_freq_beta <- ldlc.clumps_freq %>%
  inner_join(ldlc.sumstats %>% rename(ID = VARIANT) %>% select(ID, BETA), by = "ID")

# Filter by MAF >= 0.01
ldlc.instruments <- ldlc.clumps_freq_beta %>%
  filter(ALT_FREQS >= 0.01, ALT_FREQS <= 0.99)  # also exclude rare alleles > 0.99

n.ldlc <- 469680   # replace with your actual exposure GWAS sample size

ldlc.instruments <- ldlc.instruments %>%
  mutate(
    R2 = 2 * ALT_FREQS * (1 - ALT_FREQS) * BETA^2,
    F = (R2 * (n.ldlc - 2)) / (1 - R2)
  )

# Calculate summary metrics
ldlc.summary_metrics <- ldlc.instruments %>%
  summarise(
    num_snps = n(),
    cumulative_R2 = sum(R2, na.rm = TRUE),    
    mean_F = mean(F, na.rm=TRUE),
    median_F = median(F, na.rm=TRUE),
    mean_maf = mean(ALT_FREQS, na.rm = TRUE),
    mean_beta = mean(abs(BETA), na.rm = TRUE),
    overall_F = (cumulative_R2 * (n.ldlc - num_snps - 1)) / ((1 - cumulative_R2) * num_snps)
  )

kable(ldlc.summary_metrics, caption="Summary of LDL-C instruments")
```

::: {.cell-output-display}
Table: Summary of LDL-C instruments

| num_snps| cumulative_R2|   mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------:|--------:|--------:|---------:|---------:|---------:|
|      313|     0.1236359| 185.9985| 64.24422| 0.3087037| 0.0397375|   211.557|
:::
:::


  


## Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}
```
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

Matrix products: default
BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRblas.so 
LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0

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
 [1] knitr_1.49      lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
 [5] dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1    
 [9] tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] bit_4.5.0.1       gtable_0.3.6      jsonlite_1.8.9    crayon_1.5.3     
 [5] compiler_4.4.0    tidyselect_1.2.1  parallel_4.4.0    scales_1.3.0     
 [9] yaml_2.3.9        fastmap_1.2.0     R6_2.5.1          generics_0.1.3   
[13] htmlwidgets_1.6.4 munsell_0.5.1     pillar_1.10.1     tzdb_0.4.0       
[17] rlang_1.1.4       stringi_1.8.4     xfun_0.50         bit64_4.5.2      
[21] timechange_0.3.0  cli_3.6.3         withr_3.0.2       magrittr_2.0.3   
[25] digest_0.6.37     grid_4.4.0        vroom_1.6.5       rstudioapi_0.17.1
[29] hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       evaluate_1.0.1   
[33] glue_1.8.0        colorspace_2.1-1  rmarkdown_2.29    tools_4.4.0      
[37] pkgconfig_2.0.3   htmltools_0.5.8.1
```
:::
:::