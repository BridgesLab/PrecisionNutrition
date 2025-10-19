---
title: "MR Analyses of Secondary Outcomes for Total Cholesterol on Calcium Homeostasis"
author: "Dave Bridges"
date: "October 16, 2025"
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
bibliography: references.bib
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

To test if SNPs for total cholesterol GWAS identified using UK Biobank relate to other mechanistic or pathological outcomes related to calcium homeostasis and bone health.  This script can be found in /Users/davebrid/Documents/GitHub/PrecisionNutrition/Human Genetics and was most recently run on Sun Oct 19 13:05:56 2025

## Data Entry


::: {.cell}

```{.r .cell-code}
instruments.tc.file <- 'Total Cholesterol Instruments from UKBB.csv'

# loaded and renamed columns
instruments.tc <- read_csv(instruments.tc.file) |>
  rename(
    SNP                       = SP2,
    beta.exposure             = BETA,
    se.exposure               = SE,
    effect_allele.exposure    = EA,
    other_allele.exposure     = OA,
    pval.exposure             = P,
    eaf.exposure              = ALT_FREQS,
    samplesize.exposure       = N_exposure
  ) |>
  mutate(id.exposure="Total Cholesterol (UK Biobank)",
         exposure="Total Cholesterol (UK Biobank)")
```
:::


We used 370 SNPs as instruments for total cholesterol from UK Biobank.  These are found in the Total Cholesterol Instruments from UKBB.csv datafile.

## Mechanistic Outcomes

### Vitamin D Levels

This analysis is to test the hypothesis that total cholesterol impacts 25-hydroxyvitamin D levels, as cholesterol is a precursor for vitamin D synthesis in the skin.  This could positively impact calcium levels indirectly via increased vitamin D.  This was from MGI-BioVU LabWAS data [@goldsteinLabWASNovelFindings2020].


::: {.cell}

```{.r .cell-code}
gwas.vitd.file <- 'PheWeb Summary Statistics/phenocode-Vit-D.tsv.gz'
samplesize.outcome.vitd <- 12250  # sample size for MGI/BioVU for calcium from 

gwas.vitd <- read_tsv(gwas.vitd.file) |>
  mutate(ID=paste(chrom, pos, ref,alt, sep=":")) |>
  rename(
    SNP                        = ID,            # or ID if that’s the matching ID
    beta.outcome               = beta,
    se.outcome                 = sebeta,
    effect_allele.outcome      = alt,   # whichever is effect allele
    other_allele.outcome       = ref,   # whichever is other allele
    pval.outcome               = pval,
    eaf.outcome                = maf,
  ) |>
  mutate(id.outcome = "Vitamin D (MGI-BioVU LabWAS)",
         outcome = "Vitamin D (MGI-BioVU LabWAS)",
         samplesize.outcome = samplesize.outcome.vitd)  # sample size for MGI/BioVU for vitamin D)
library(TwoSampleMR)

vitd.data <- harmonise_data(instruments.tc, gwas.vitd, action = 2)
vitd.data_steiger <- steiger_filtering(vitd.data)
vitd.mr <- mr(vitd.data_steiger,
                         method_list = c("mr_ivw", 
                                         "mr_egger_regression",
                                         "mr_weighted_median", 
                                         "mr_weighted_mode"))

vitd.mr |> select(-starts_with('id')) |> 
  kable(caption="MR Results for Total Cholesterol - Vitamin D Analysis",
        digits=c(0,0,0,0,3,3,99))
```

::: {.cell-output-display}


Table: MR Results for Total Cholesterol - Vitamin D Analysis

|outcome                      |exposure                       |method                    | nsnp|      b|    se|       pval|
|:----------------------------|:------------------------------|:-------------------------|----:|------:|-----:|----------:|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Inverse variance weighted |  280| -0.063| 0.033| 0.05362791|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |MR Egger                  |  280| -0.070| 0.053| 0.19017784|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted median           |  280| -0.005| 0.052| 0.92221486|
|Vitamin D (MGI-BioVU LabWAS) |Total Cholesterol (UK Biobank) |Weighted mode             |  280| -0.012| 0.051| 0.82052248|


:::

```{.r .cell-code}
ggplot(vitd.mr, aes(y=method,x=b)) +
  geom_point() +
  geom_errorbar(aes(xmin=b-1.96*se, xmax=b+1.96*se), width=0.2) +
  theme_classic(base_size=16) +
  labs(title="25-Hydroxyvitamin D Levels (LabWAS)",
       y="",
       x="Effect Size (Beta)") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") 
```

::: {.cell-output-display}
![](figures/mr-vitd-1.png){width=672}
:::
:::


## Pathological Outcomes

### Bone Mineral Density

#### GEFOS Consortium Lumbar Spine BMD

Lumbar spine tends to reflect trabecular bone, which is metabolically active, hormonally responsive, and more sensitive to lipid and endocrine perturbations. Multiple epidemiologic studies report inverse associations between total cholesterol and LS-BMD, particularly among older adults and postmenopausal women [@huAssociationsSerumTotal2023;@fangNegativeAssociationTotal2022].  We used the pooled (not sex-specific) estimates from the GEFOS consortium's pooled meta-analysis [@estradaGenomewideMetaanalysisIdentifies2012].



::: {.cell}

```{.r .cell-code}
gwas.bmd.file <- 'PheWeb Summary Statistics/GEFOS2_LSBMD_POOLED_GC.txt.gz'
samplesize.outcome.bmd <- 32961  # sample size from estradaGenomewideMetaanalysisIdentifies2012

gwas.bmd <- read_tsv(gwas.bmd.file) |>
  rename(
    beta.outcome               = Effect,
    se.outcome                 = StdErr,
    pval.outcome               = `P-value`,
    eaf.outcome                = Freq1,
  ) |>
  mutate(id.outcome = "bmd",
         outcome = "LS-BMD (GEFOS)",
         effect_allele.outcome = toupper(Allele1),   # whichever is effect allele
         other_allele.outcome = toupper(Allele2),   # whichever is other allele
         samplesize.outcome = samplesize.outcome.bmd)  # sample size for MGI/BioVU for vitamin D)

#runs very slow so commented out after saving the SNP file 
library(biomaRt)
#snp.mart <- useEnsembl(biomart="snp", dataset="hsapiens_snp")
#snp.data <- getBM(attributes=c("refsnp_id", "chr_name", "chrom_start", "allele"),
#                  filters="snp_filter",
#                  values=gwas.bmd$MarkerName,
#                  mart=snp.mart)
gefos.snp.datafile <- "GEFOS-2012 SNP ID File.csv"
#write_csv(snp.data, gefos.snp.datafile)
bmd.snp.data <- read_csv(gefos.snp.datafile) |>
  mutate(SNP=paste(chr_name, chrom_start,
                        substr(allele, 1, 1), sep=":")) #just picked hte first alt allele

gwas.bmd.combined <- 
  left_join(gwas.bmd,bmd.snp.data, by=c("MarkerName"="refsnp_id")) 
  
#no harmonised SNPs
bmd.data <- harmonise_data(instruments.tc, gwas.bmd.combined, action = 2)
bmd.data_steiger <- steiger_filtering(bmd.data)
bmd.mr <- mr(bmd.data_steiger,
                         method_list = c("mr_ivw", 
                                         "mr_egger_regression",
                                         "mr_weighted_median", 
                                         "mr_weighted_mode"))

bmd.mr |> dplyr::select(-starts_with('id')) |> 
  kable(caption="MR Results for Total Cholesterol - Lumbar Spine BMD (GEFOS - 2012)",
        digits=c(0,0,0,0,3,3,99))
```

::: {.cell-output-display}


Table: MR Results for Total Cholesterol - Lumbar Spine BMD (GEFOS - 2012)

|SNP |MarkerName |Allele1 |Allele2 | eaf.outcome| FreqSE| beta.outcome| se.outcome| pval.outcome|Direction |outcome |effect_allele.outcome |other_allele.outcome | samplesize.outcome|chr_name | chrom_start|allele | CHR| POS|effect_allele.exposure |other_allele.exposure | beta.exposure| se.exposure| pval.exposure| eaf.exposure| samplesize.exposure| R2|  F|exposure |
|:---|:----------|:-------|:-------|-----------:|------:|------------:|----------:|------------:|:---------|:-------|:---------------------|:--------------------|------------------:|:--------|-----------:|:------|---:|---:|:----------------------|:---------------------|-------------:|-----------:|-------------:|------------:|-------------------:|--:|--:|:--------|


:::

```{.r .cell-code}
# ggplot(bmd.mr, aes(y=method,x=b)) +
#   geom_point() +
#   geom_errorbar(aes(xmin=b-1.96*se, xmax=b+1.96*se), width=0.2) +
#   theme_classic(base_size=16) +
#   labs(title="Bone Minderal Density (GEFOS)",
#        y="",
#        x="Effect Size (Beta)") +
#   geom_vline(xintercept=0, linetype="dashed", color = "red") 
```
:::

This second analysis is from lumbar spine data from the GEFOS 2015 release ([@zhengWholegenomeSequencingIdentifies2015])


::: {.cell}

```{.r .cell-code}
gwas.bmd.file.2015 <- 'PheWeb Summary Statistics/ls2stu.MAF0_.005.pos_.out_.gz'
samplesize.outcome.bmd.2015 <- 53236  # sample size from zhengWholegenomeSequencingIdentifies2015

gwas.bmd.2015 <- read_tsv(gwas.bmd.file.2015) |>
  rename(
    beta.outcome               = beta,
    se.outcome                 = se,
    effect_allele.outcome      = reference_allele,   # whichever is effect allele
    other_allele.outcome       = other_allele,   # whichever is other allele
    pval.outcome               = `p-value`,
    eaf.outcome                = eaf,
  ) |>
  mutate(SNP=paste(chromosome,position,effect_allele.outcome,other_allele.outcome,sep=":"),
         id.outcome = "bmd",
         outcome = "LS-BMD (GEFOS)",
         samplesize.outcome = samplesize.outcome.bmd.2015)  # sample size for MGI/BioVU for vitamin D)

bmd.data.2015 <- harmonise_data(instruments.tc, gwas.bmd.2015, action = 2)
bmd.data_steiger.2015 <- steiger_filtering(bmd.data.2015)
bmd.mr.2015 <- mr(bmd.data_steiger.2015,
                         method_list = c("mr_ivw", 
                                         "mr_egger_regression",
                                         "mr_weighted_median", 
                                         "mr_weighted_mode"))

bmd.mr.2015 |> 
  dplyr::select(-starts_with('id')) |> 
  kable(caption="MR Results for Total Cholesterol - Lumbar Spine BMD (GEFOS - 2015)",
        digits=c(0,0,0,0,3,3,99))
```

::: {.cell-output-display}


Table: MR Results for Total Cholesterol - Lumbar Spine BMD (GEFOS - 2015)

|outcome        |exposure                       |method                    | nsnp|      b|    se|      pval|
|:--------------|:------------------------------|:-------------------------|----:|------:|-----:|---------:|
|LS-BMD (GEFOS) |Total Cholesterol (UK Biobank) |Inverse variance weighted |   84| -0.061| 0.045| 0.1744722|
|LS-BMD (GEFOS) |Total Cholesterol (UK Biobank) |MR Egger                  |   84| -0.060| 0.071| 0.4029161|
|LS-BMD (GEFOS) |Total Cholesterol (UK Biobank) |Weighted median           |   84|  0.042| 0.063| 0.5056823|
|LS-BMD (GEFOS) |Total Cholesterol (UK Biobank) |Weighted mode             |   84|  0.009| 0.061| 0.8871697|


:::

```{.r .cell-code}
ggplot(bmd.mr.2015, aes(y=method,x=b)) +
  geom_point() +
  geom_errorbar(aes(xmin=b-1.96*se, xmax=b+1.96*se), width=0.2) +
  theme_classic(base_size=16) +
  labs(title="Bone Mineral Density (GEFOS)",
       y="",
       x="Effect Size (Beta)") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") 
```

::: {.cell-output-display}
![](figures/mr-bmd-2015-1.png){width=672}
:::
:::


### Fracture Risk

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
 [1] biomaRt_2.64.0     TwoSampleMR_0.6.22 knitr_1.50         lubridate_1.9.4   
 [5] forcats_1.0.1      stringr_1.5.2      dplyr_1.1.4        purrr_1.1.0       
 [9] readr_2.1.5        tidyr_1.3.1        tibble_3.3.0       ggplot2_4.0.0     
[13] tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] KEGGREST_1.48.1         gtable_0.3.6            httr2_1.2.1            
 [4] xfun_0.53               htmlwidgets_1.6.4       psych_2.5.6            
 [7] Biobase_2.68.0          lattice_0.22-7          tzdb_0.5.0             
[10] vctrs_0.6.5             tools_4.5.1             generics_0.1.4         
[13] curl_7.0.0              stats4_4.5.1            parallel_4.5.1         
[16] AnnotationDbi_1.70.0    RSQLite_2.4.3           blob_1.2.4             
[19] pkgconfig_2.0.3         data.table_1.17.8       dbplyr_2.5.1           
[22] RColorBrewer_1.1-3      S7_0.2.0                S4Vectors_0.46.0       
[25] GenomeInfoDbData_1.2.14 lifecycle_1.0.4         compiler_4.5.1         
[28] farver_2.1.2            progress_1.2.3          Biostrings_2.76.0      
[31] mnormt_2.1.1            GenomeInfoDb_1.44.3     htmltools_0.5.8.1      
[34] yaml_2.3.10             pillar_1.11.1           crayon_1.5.3           
[37] cachem_1.1.0            nlme_3.1-168            tidyselect_1.2.1       
[40] digest_0.6.37           stringi_1.8.7           labeling_0.4.3         
[43] fastmap_1.2.0           grid_4.5.1              cli_3.6.5              
[46] magrittr_2.0.4          withr_3.0.2             filelock_1.0.3         
[49] rappdirs_0.3.3          prettyunits_1.2.0       UCSC.utils_1.4.0       
[52] scales_1.4.0            bit64_4.6.0-1           timechange_0.3.0       
[55] XVector_0.48.0          httr_1.4.7              rmarkdown_2.30         
[58] bit_4.6.0               png_0.1-8               hms_1.1.4              
[61] memoise_2.0.1           evaluate_1.0.5          IRanges_2.42.0         
[64] BiocFileCache_2.16.2    rlang_1.1.6             Rcpp_1.1.0             
[67] glue_1.8.0              DBI_1.2.3               xml2_1.4.0             
[70] BiocGenerics_0.54.1     rstudioapi_0.17.1       vroom_1.6.6            
[73] jsonlite_2.0.0          R6_2.6.1                plyr_1.8.9             
```


:::
:::

