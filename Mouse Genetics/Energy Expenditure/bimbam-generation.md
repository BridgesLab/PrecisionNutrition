---
title: "Generation of BIMBAM files for GEMMA Analyses"
author: "Dave Bridges"
date: "February 11, 2022"
output:
  html_document:
    highlight: tango
    keep_md: yes
    number_sections: yes
    toc: yes
  pdf_document:
    highlight: tango
    keep_tex: yes
    number_sections: yes
    toc: yes
---



# Purpose

Use this file to generate bimbam files for GEMMA analyses

# Raw Data

Downloaded SNPs from http://gn1.genenetwork.org/dbdoc/BXDGeno.html


```r
bxd.genotypes <- 'BXD-geno_8-9-17.txt'
library(readr)
bxd.genotype.data <- read_tsv(bxd.genotypes, skip=20,
                              col_types = cols(
  .default = col_character(),
  Chr = col_factor(levels=NULL),
  cM = col_double(),
  Mb = col_double()
)) %>%  
  mutate(`5XFAD`='B') %>% #strain 5x FAD is set to homozygous B allele
  mutate(`D2Gpnmb`='D') %>% #strain 5x FAD is set to homozygous B allele
  mutate(`D2`='D')  #strain 5x FAD is set to homozygous B allele

strain.labels <- 'bxd-strain-labels.txt'
colnames(bxd.genotype.data) %>% tail(-4) %>% write(strain.labels)
```

These data can be found in **/Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/Other Published Data/Systems Biology**.  This script was most recently updated on **Fri Feb 11 12:56:39 2022**.

This analysis uses the BXD genotyping file at BXD-geno_8-9-17.txt which had 7320 markers for 205 strains of mice.

## SNP Annotation File

The SNP annotation file must be in the order SNP - bp - chromosome


```r
snp.annotation.file <- 'snp-annotation.bimbam'
bxd.genotype.data %>%
  mutate(Base=Mb*1000000) %>% #convert Mb to actual base position
  select(Locus,Base,Chr) %>%
  write_csv(snp.annotation.file, col_names = FALSE)
```

## Genotype File

The genotype file must be in the order SNP - Major Allele - Minor Allele - Genotype1 ... GenotypeX.  This genotype file contains all genotypes, including those with no data.


```r
genotype.file <- 'bxd-genotype.bimbam' 
bxd.genotype.data %>%
  mutate(Major = 'B',
         Minor = 'D') %>%
  select(Locus,Major,Minor,
         starts_with('BXD'),
         starts_with('5XFAD'),
         starts_with('D2')) %>% 
  write_csv(genotype.file, col_names=FALSE)
```

# Session Information


```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS  10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] readr_2.1.1 dplyr_1.0.7 tidyr_1.1.4 knitr_1.37 
## 
## loaded via a namespace (and not attached):
##  [1] pillar_1.6.4     bslib_0.3.1      compiler_4.0.2   jquerylib_0.1.4 
##  [5] tools_4.0.2      bit_4.0.4        digest_0.6.29    jsonlite_1.7.2  
##  [9] evaluate_0.14    lifecycle_1.0.1  tibble_3.1.6     pkgconfig_2.0.3 
## [13] rlang_0.4.12     DBI_1.1.2        parallel_4.0.2   yaml_2.2.1      
## [17] xfun_0.29        fastmap_1.1.0    stringr_1.4.0    generics_0.1.1  
## [21] vctrs_0.3.8      sass_0.4.0       hms_1.1.1        bit64_4.0.5     
## [25] tidyselect_1.1.1 glue_1.6.0       R6_2.5.1         fansi_1.0.0     
## [29] vroom_1.5.7      rmarkdown_2.11   purrr_0.3.4      tzdb_0.2.0      
## [33] magrittr_2.0.1   ellipsis_0.3.2   htmltools_0.5.2  assertthat_0.2.1
## [37] utf8_1.2.2       stringi_1.7.6    crayon_1.4.2
```
