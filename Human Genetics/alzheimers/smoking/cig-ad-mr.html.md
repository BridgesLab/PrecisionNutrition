---
title: "MR Analyses of Cigarettes on AD Risk"
author: "Dave Bridges and Katie Kittell"
date: "January 19, 2026"
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

To test if SNPs for sleep duration identified using Dashti *et al.* relate to AD risk.  This script can be found in /Volumes/BridgesLab/Kittell/MR - Sleep:AD and was most recently run on Tue May 19 10:37:27 2026

## Data Entry




::: {.cell}

```{.r .cell-code}
library(TwoSampleMR)

# Cigarettes
instruments.file <- '/Volumes/BridgesLab/Kittell/MR - Sleep:AD/Instruments - cig.liu - Sleep.csv'
instruments <- read_csv(instruments.file) |>
  mutate(samplesize.exposure=249752) |>
  mutate(exposure = "Cigarettes")


#OpenGWAS ID for Bellenguez et al. (2022)
ad_outcome_id <- "ebi-a-GCST90027158"

#Extract the outcome data for your specific SNPs
ad.outcome <- extract_outcome_data(
    snps = instruments$SNP,      # Your list of instrument SNPs
    outcomes = ad_outcome_id,     # The dataset ID
    proxies = TRUE,               # (Optional) Look for proxies if a SNP is missing
    rsq = 0.8)                     # (Optional) LD threshold for proxies


samplesize.ad.outcome <- 487511 
```
:::





We used 23 SNPs as instruments for Cigarettes  These are found in the /Volumes/BridgesLab/Kittell/MR - Sleep:AD/Instruments - cig.liu - Sleep.csv datafile.





::: {.cell}

```{.r .cell-code}
ad.data <- harmonise_data(instruments, ad.outcome, action = 2)
ad.data.steiger <- steiger_filtering(ad.data) #need to do steiger filtering

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# Create GRanges object from your SNP data
snps <- GRanges(seqnames = paste0("chr", ad.data.steiger$chr),
                ranges = IRanges(start = as.numeric(ad.data.steiger$pos),
                                end = as.numeric(ad.data.steiger$pos)),
                rsid = ad.data.steiger$SNP)

# Get all genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb)

# Find nearest gene for EVERY SNP
nearest_idx <- nearest(snps, genes)
distances <- distance(snps, genes[nearest_idx])

# Get gene symbols
gene_ids <- genes[nearest_idx]$gene_id
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID")

# Combine results
snp_annotation_results <- data.frame(
  snp = snps$rsid,
  chr = seqnames(snps),
  pos = start(snps),
  nearest_gene = gene_symbols,
  distance = distances
)

#instrument strength
ad.data.steiger.annot <- 
    ad.data.steiger %>%
  left_join(snp_annotation_results |> dplyr::select(snp, nearest_gene, distance),
            by = c("SNP"="snp")) |>
  filter(mr_keep==TRUE) |>
  mutate(
    R2.exposure = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
    F.exposure = (R2.exposure * (samplesize.exposure - 2)) / (1 - R2.exposure)
  )

ad.data.steiger.annot.summary <- ad.data.steiger.annot %>%
  summarise(
    num_snps = n(),
    samplesize.exposure = dplyr::first(samplesize.exposure),
    cumulative_R2 = sum(R2.exposure, na.rm = TRUE),
    mean_F = mean(F.exposure, na.rm = TRUE),
    median_F = median(F.exposure, na.rm = TRUE),
    mean_maf = mean(eaf.exposure, na.rm = TRUE),
    mean_beta = mean(abs(beta.exposure), na.rm = TRUE)
  ) |>
  mutate(overall_F = (cumulative_R2 * (samplesize.exposure - num_snps - 1)) / 
                     ((1 - cumulative_R2) * num_snps))

library(knitr)
kable(ad.data.steiger.annot.summary, caption="Summary of Cigarettes instruments after harmonisation for AD analysis")
```

::: {.cell-output-display}


Table: Summary of Cigarettes instruments after harmonisation for AD analysis

| num_snps| samplesize.exposure| cumulative_R2|   mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------------:|-------------:|--------:|--------:|---------:|---------:|---------:|
|       22|              249752|     0.0348919| 399.0924| 153.8581| 0.3728955| 0.0582883|  410.3887|


:::

```{.r .cell-code}
# write to CSV
write.csv(ad.data.steiger.annot.summary,
          file = file.path("exposure_tables", "exposure_table_Cigarettes.csv"),
          row.names = FALSE)
```
:::




After harmonization, we were left with 23 SNPs. 1 SNPs removed due to allelle ambiguity.  After Steiger filtering we had a total number of 22 SNPs for analysis.  

## MR Analysis

We analyzed this potential association with a series of MR methods.




::: {.cell}

```{.r .cell-code}
ad.mr <- mr(ad.data.steiger,
                         method_list = c("mr_ivw_mre",
                                         "mr_ivw_fe",
                                         "mr_raps",
                                         "mr_egger_regression",
                                         "mr_weighted_median", 
                                         "mr_weighted_mode"))

mr_pleiotropy_test(ad.data.steiger) |>
  dplyr::select(-starts_with('id')) |> 
  kable(caption="MR Pleiotropy Results for Cigarettes-AD Analysis")
```

::: {.cell-output-display}


Table: MR Pleiotropy Results for Cigarettes-AD Analysis

|outcome                                                |exposure   | egger_intercept|        se|      pval|
|:------------------------------------------------------|:----------|---------------:|---------:|---------:|
|Alzheimer's disease &#124;&#124; id:ebi-a-GCST90027158 |Cigarettes |       0.0114142| 0.0067915| 0.1083789|


:::

```{.r .cell-code}
mr_heterogeneity(ad.data.steiger) |>
  dplyr::select(-starts_with('id')) |> 
    mutate(
    I2 = pmax(0, (Q - Q_df) / Q) * 100 # 
  ) |>
  kable(caption="MR Heterogeneity Results for Cigarettes-AD Analysis",
        digits=c(0,0,0,3,3,99))
```

::: {.cell-output-display}


Table: MR Heterogeneity Results for Cigarettes-AD Analysis

|outcome                                                |exposure   |method                    |      Q| Q_df|       Q_pval| I2|
|:------------------------------------------------------|:----------|:-------------------------|------:|----:|------------:|--:|
|Alzheimer's disease &#124;&#124; id:ebi-a-GCST90027158 |Cigarettes |MR Egger                  | 72.926|   20| 6.020378e-08| 73|
|Alzheimer's disease &#124;&#124; id:ebi-a-GCST90027158 |Cigarettes |Inverse variance weighted | 83.225|   21| 2.317303e-09| 75|


:::

```{.r .cell-code}
#M-PRESSO has to be run separately
library(MRPRESSO)

ad.mr_presso_results <- mr_presso(
  BetaOutcome = "beta.outcome",      # Column name for outcome betas
  BetaExposure = "beta.exposure",    # Column name for exposure betas
  SdOutcome = "se.outcome",          # Column name for outcome SEs
  SdExposure = "se.exposure",        # Column name for exposure SEs
  data = ad.data.steiger,                  # Your dataset
  NbDistribution = 2000,              # Number of distributions (default 1000)
  SignifThreshold = 0.05,             # Significance threshold
  OUTLIERtest = TRUE,                 # Perform outlier test
  DISTORTIONtest = TRUE,              # Perform distortion test
)

library(forcats)
ad.mr_presso_results$`Main MR results` |>
  dplyr::select(`MR Analysis`, `Causal Estimate`, Sd, `P-value`) |>
  dplyr::rename(method = `MR Analysis`,
         b = `Causal Estimate`,
         se = Sd,
         pval = `P-value`) |>
  mutate(method = fct_recode(method,
                             "MR-PRESSO (Raw)"="Raw",
                             "MR-PRESSO (Outlier-corrected)"="Outlier-corrected")) -> ad.mr_presso_results_df

ad.mr.mrpresso <- bind_rows(as_tibble(ad.mr), as_tibble(ad.mr_presso_results_df))|>
  fill(id.exposure, id.outcome, exposure, outcome,nsnp,.direction="down")

Cigarettes_results_df <- ad.mr.mrpresso |> 
  dplyr::select(-starts_with('id')) |> 
  kable(caption="MR Results for Cigarettes-AD Analysis",
        digits=c(0,0,0,0,3,3,99))

# write csv
write.csv(Cigarettes_results_df, 
          file = file.path("AD_MR_results","ad_Cigarettes_results.csv"), 
                           row.names = FALSE)

ggplot(ad.mr.mrpresso, aes(y=method,x=b)) +
  geom_point() +
  geom_errorbar(aes(xmin=b-1.96*se, xmax=b+1.96*se), width=0.2) +
  theme_classic(base_size=16) +
  theme(axis.text.y = element_text(size = 10)) +
  labs(title="AD Risk (Cigarettes)",
       y="",
       x="Effect Size (Beta)") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") 
```

::: {.cell-output-display}
![](figures/ad-mr-Cigarettes-1.png){width=672}
:::

```{.r .cell-code}
beta.ivw <- filter(ad.mr, method == "Inverse variance weighted (multiplicative random effects)")$b
se.ivw <- filter(ad.mr, method == "Inverse variance weighted (multiplicative random effects)")$se
p.ivw <- filter(ad.mr, method == "Inverse variance weighted (multiplicative random effects)")$pval

# =============================================================================
# Bayesian Analysis of MR Results
# =============================================================================

# Observational hypothesis
obs_OR <- 0.8
obs_mean <- log(obs_OR)  # -0.223
obs_sd <- 0.1  # uncertainty on log OR scale

# =============================================================================
# BAYES FACTORS
# =============================================================================

# BF1: Observational effect vs Null
# ----------------------------------
# H_obs: β ~ N(log(0.8), 0.1^2)
# H_null: β = 0

likelihood_obs <- dnorm(beta.ivw, 
                        mean = obs_mean, 
                        sd = sqrt(se.ivw^2 + obs_sd^2))

likelihood_null <- dnorm(beta.ivw, 
                         mean = 0, 
                         sd = se.ivw)

BF_obs_vs_null <- likelihood_obs / likelihood_null

# BF2: Any effect vs Null (using standard Wakefield prior)
# ---------------------------------------------------------
# H_any: β ~ N(0, 0.15^2) - allows effect in either direction
# H_null: β = 0

W <- 0.15^2  # standard prior variance for epidemiological effects
V <- se.ivw^2

BF_any_vs_null <- sqrt((W + V) / V) * 
                  exp(-(beta.ivw^2 / 2) * (W / (V * (W + V))))

# =============================================================================
# POSTERIOR PROBABILITIES
# =============================================================================

# Use weakly informative prior centered at 0
prior_mean <- 0
prior_sd <- 0.2  # weak prior allowing broad range of effects

# Posterior distribution
post_precision <- 1/se.ivw^2 + 1/prior_sd^2
post_var <- 1/post_precision
post_mean <- post_var * (beta.ivw/se.ivw^2 + prior_mean/prior_sd^2)
post_sd <- sqrt(post_var)

# Posterior probability 1: Effect is protective (β < 0)
p_protective <- pnorm(0, mean = post_mean, sd = post_sd, lower.tail = TRUE)

# Posterior probability 2: Effect within observational range
# Define range as OR 0.8 ± 0.1 on log scale
range_lower <- obs_mean - obs_sd  # log(0.8) - 0.1 = -0.323
range_upper <- obs_mean + obs_sd  # log(0.8) + 0.1 = -0.123

p_in_obs_range <- pnorm(range_upper, mean = post_mean, sd = post_sd) - 
                  pnorm(range_lower, mean = post_mean, sd = post_sd)


# Convert posterior to OR scale for interpretation
post_OR <- exp(post_mean)
post_OR_lower <- exp(post_mean - 1.96*post_sd)
post_OR_upper <- exp(post_mean + 1.96*post_sd)
```
:::




### Scatter Plot

The scatter plot visualizes the individual SNPs' effects on the exposure (Cigarettes) and outcome (AD), along with the overall MR estimate as a regression line.  Nearest genes are included for annotation.




::: {.cell}

```{.r .cell-code}
library(ggrepel)
ggplot(ad.data.steiger.annot, aes(x=beta.exposure, y=beta.outcome)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome,
                    ymax = beta.outcome + 1.96*se.outcome),
                alpha=0.5) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96*se.exposure,
                    xmax = beta.exposure + 1.96*se.exposure),
                alpha=0.5) +
  geom_smooth(method="lm",se=F) +
  geom_text_repel(aes(label=nearest_gene),
                  max.overlaps = 15) + 
  theme_classic(base_size=16) +
  labs(x="Exposure Estimate (Cigarettes)", 
       y="Outcome Estimate (AD)", 
       title="") 
```

::: {.cell-output-display}
![](figures/mr-scatter-Cigarettes-1.png){width=672}
:::
:::




### Funnel Plot

Funnel plots help visualize heterogeneity and potential directional horizontal pleiotropy




::: {.cell}

```{.r .cell-code}
single_snp_results <- mr_singlesnp(ad.data.steiger.annot)

# Get overall IVW estimate for the vertical line
ivw_beta <- ad.mr |> filter(method=="Inverse variance weighted (multiplicative random effects)") |> pull(b)

# Determine y-range based on your data
y_min <- 0
y_max <- max(single_snp_results$se^{-1}) * 1.1  # 10% padding above max precision

# Generate a fine grid of precision values
precision_grid <- seq(y_min, y_max, length.out = 1000)

# Compute boundaries: ivw_beta ± 1.96 / precision
lower_bound <- ivw_beta - 1.96 / precision_grid
upper_bound <- ivw_beta + 1.96 / precision_grid

# Create data frame for boundaries
bounds_df <- data.frame(precision = precision_grid, lower = lower_bound, upper = upper_bound)

# Plot
ggplot(single_snp_results, aes(x = b, y = 1/se)) +
  # Scatter points for each SNP
  geom_point(size = 1) +
  # Vertical line at IVW estimate
  geom_vline(xintercept = ivw_beta, linetype = "solid", color = "#ff7f0e", size = 1) +
  # Curved pseudo-95% CI boundaries (the cone)
  geom_line(data = bounds_df, aes(x = lower, y = precision), linetype = "dashed") +
  geom_line(data = bounds_df, aes(x = upper, y = precision), linetype = "dashed") +
  # Customize axes and labels
  labs(
    x = "Estimate (Beta-IVW)",
    y = "Precision (1/Standard Error)",
    title = ""
  ) +
  # Apply clean theme and limit y to >=0
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, y_max), xlim = c(min(single_snp_results$b), max(single_snp_results$b)))  # Adjust x-limits for visibility
```

::: {.cell-output-display}
![](figures/Cigarettes-ad-funnel-plot-1.png){width=672}
:::
:::




### Leave-one-out Analysis

This analysis shows estimates as we sequentially remove one SNP at a time, to determine if the association is strongly influenced by a single SNP




::: {.cell}

```{.r .cell-code}
# LOO using IVW
loo_res <- mr_leaveoneout(ad.data.steiger.annot)
loo_res |> 
  mutate(diff = b - filter(ad.mr, method=="Inverse variance weighted (multiplicative random effects)")$b) |>
  arrange(-abs(diff)) |>
  head() |>
  dplyr::select(SNP,diff,b,se,p) |>
  kable(caption="Leave-One-Out Results (IVW method) for Influential SNPs",
        digits=c(0,5,5,5,5))
```

::: {.cell-output-display}


Table: Leave-One-Out Results (IVW method) for Influential SNPs

|SNP        |     diff|        b|      se|       p|
|:----------|--------:|--------:|-------:|-------:|
|rs11852372 |  0.06242| -0.00242| 0.08292| 0.97671|
|rs73229090 | -0.02375| -0.08859| 0.05093| 0.08195|
|rs632811   | -0.02331| -0.08815| 0.05002| 0.07803|
|rs56113850 |  0.00986| -0.05498| 0.06832| 0.42098|
|rs2273500  |  0.00885| -0.05600| 0.06385| 0.38048|
|rs4785587  | -0.00707| -0.07191| 0.06284| 0.25246|


:::

```{.r .cell-code}
ggplot(loo_res, aes(x = reorder(SNP, -b), y = b)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin = b - 1.96*se, ymax = b + 1.96*se), width = 0.01 ,alpha=0.5) +
  coord_flip() +
  labs(x = "SNP Removed", y = "Estimate (Beta-IVW; leave-one-out)") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic(base_size=16) +
  theme(axis.text.y = element_text(size = 4)) 
```

::: {.cell-output-display}
![](figures/Cigarettes-ad-mr-loo-1.png){width=672}
:::

```{.r .cell-code}
#
```
:::







## Interpretation

"Our analysis yielded a non-significant estimate for the effect of sleep duration on AD risk ($\beta = -0.0648418, p =0.2958587$). However, the direction of the effect was consistent with our biological hypothesis (where increased sleep duration reduces risk). Given the low heritability of the sleep duration trait ($R^2=0.0348919$) and the resulting low precision ($SE=0.0620284$), we performed a Bayesian analysis to assess the strength of this null result. 

## Bayesian Analysis of Mendelian Randomization Results

We conducted Bayesian analyses to evaluate the evidence for a causal protective effect of sleep duration on Alzheimer's disease risk.

### Bayes Factors

We computed Bayes factors comparing two alternative hypotheses against the null hypothesis of no effect:

**Observational effect hypothesis** $\beta \sim N(log(0.8), 0.1^2)$ vs null: BF = 0.37. This indicates that the data favor the null hypothesis by approximately 2.7:1 odds.

**Any effect hypothesis** (allowing effects in either direction) vs null: 
BF = 1.64, suggesting weak evidence for an effect.

### Posterior Probabilities

Using a weakly informative prior $\beta \sim N(log(0), 0.2^2)$, we estimated the posterior distribution of the causal effect. The posterior mean effect was $\beta$ = -0.059 (95% CI: -0.175 to 0.057), corresponding to an odds ratio of 0.94 (95% CI: 0.84 to 1.06).

The posterior probability that sleep duration has a protective effect ($\beta < 0$) was 84.1%. The posterior probability that the true effect lies within the range suggested by observational studies (OR $0.8 \pm 0.1$) was 14%.

### Interpretation

While our MR estimates are consistent with a modest protective effect, the evidence is weak and does not strongly support the observational findings.


## References

::: {#refs}
:::

## Session Information




::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.4.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Detroit
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggrepel_0.9.8                          
 [2] MRPRESSO_1.0                           
 [3] org.Hs.eg.db_3.20.0                    
 [4] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
 [5] GenomicFeatures_1.58.0                 
 [6] AnnotationDbi_1.68.0                   
 [7] Biobase_2.66.0                         
 [8] GenomicRanges_1.58.0                   
 [9] GenomeInfoDb_1.42.3                    
[10] IRanges_2.40.1                         
[11] S4Vectors_0.44.0                       
[12] BiocGenerics_0.52.0                    
[13] TwoSampleMR_0.6.29                     
[14] knitr_1.51                             
[15] lubridate_1.9.5                        
[16] forcats_1.0.1                          
[17] stringr_1.6.0                          
[18] dplyr_1.2.0                            
[19] purrr_1.2.1                            
[20] readr_2.2.0                            
[21] tidyr_1.3.2                            
[22] tibble_3.3.1                           
[23] ggplot2_4.0.2                          
[24] tidyverse_2.0.0                        

loaded via a namespace (and not attached):
 [1] DBI_1.3.0                   mnormt_2.1.2               
 [3] bitops_1.0-9                gridExtra_2.3              
 [5] rlang_1.1.7                 magrittr_2.0.4             
 [7] otel_0.2.0                  matrixStats_1.5.0          
 [9] compiler_4.4.2              RSQLite_2.4.6              
[11] mgcv_1.9-4                  png_0.1-9                  
[13] vctrs_0.7.1                 httpcode_0.3.0             
[15] pkgconfig_2.0.3             crayon_1.5.3               
[17] fastmap_1.2.0               XVector_0.46.0             
[19] labeling_0.4.3              Rsamtools_2.22.0           
[21] rmarkdown_2.30              tzdb_0.5.0                 
[23] UCSC.utils_1.2.0            bit_4.6.0                  
[25] xfun_0.56                   zlibbioc_1.52.0            
[27] cachem_1.1.0                jsonlite_2.0.0             
[29] blob_1.3.0                  DelayedArray_0.32.0        
[31] mr.raps_0.4.3               BiocParallel_1.40.2        
[33] psych_2.6.1                 parallel_4.4.2             
[35] R6_2.6.1                    stringi_1.8.7              
[37] RColorBrewer_1.1-3          rtracklayer_1.66.0         
[39] Rcpp_1.1.1                  SummarizedExperiment_1.36.0
[41] splines_4.4.2               Matrix_1.7-4               
[43] timechange_0.4.0            tidyselect_1.2.1           
[45] rstudioapi_0.18.0           abind_1.4-8                
[47] yaml_2.3.12                 codetools_0.2-20           
[49] curl_7.0.0                  lattice_0.22-9             
[51] plyr_1.8.9                  withr_3.0.2                
[53] KEGGREST_1.46.0             S7_0.2.1                   
[55] evaluate_1.0.5              Biostrings_2.74.1          
[57] pillar_1.11.1               MatrixGenerics_1.18.1      
[59] nortest_1.0-4               ieugwasr_1.1.0.9000        
[61] generics_0.1.4              vroom_1.7.0                
[63] RCurl_1.98-1.17             hms_1.1.4                  
[65] rootSolve_1.8.2.4           scales_1.4.0               
[67] glue_1.8.0                  tools_4.4.2                
[69] BiocIO_1.16.0               data.table_1.18.2.1        
[71] rsnps_0.6.1                 GenomicAlignments_1.42.0   
[73] XML_3.99-0.22               grid_4.4.2                 
[75] nlme_3.1-168                GenomeInfoDbData_1.2.13    
[77] restfulr_0.0.16             cli_3.6.5                  
[79] S4Arrays_1.6.0              gtable_0.3.6               
[81] digest_0.6.39               crul_1.6.0                 
[83] SparseArray_1.6.2           rjson_0.2.23               
[85] htmlwidgets_1.6.4           farver_2.1.2               
[87] memoise_2.0.1               htmltools_0.5.9            
[89] lifecycle_1.0.5             httr_1.4.8                 
[91] bit64_4.6.0-1              
```


:::
:::
