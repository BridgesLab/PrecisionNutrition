---
title: "Analysis of Cholesterol QTL Region on C"
author: "Dave Bridges"
format: html
execute:
  keep-md: true
---





# Data Entry

Downloaded the genomic region from Genome MUSTER from 162,000,000 to 180,000,000 for all the founder strains.  This script can be found in /Users/davebrid/Documents/GitHub/PrecisionNutrition/Mouse Genetics/Cholesterol/GWAS and was most recently run on Sun Oct 15 16:27:34 2023.


::: {.cell}

:::


# Analysis

Based on the effect estimation the allele is present in the 129S1/SvImJ but not the other strains, particularly A/J.  We ran GenomeMuster to only provide variants that are different between these strains.

# Annotation

The MVAR annotation didnt come out from genome Muster so used https://genome.ucsc.edu/cgi-bin/hgVai on the known rs id's from these strains.


::: {.cell}

:::


## Analysis of Annotation


::: {.cell}
::: {.cell-output-display}
Table: Variants by annotation type

|Annotation                         |      n|
|:----------------------------------|------:|
|intron_variant                     | 123670|
|upstream_gene_variant              |  29873|
|downstream_gene_variant            |  28047|
|intergenic_variant                 |  27345|
|NMD_transcript_variant             |   6724|
|non_coding_transcript_exon_variant |   2606|
|3_prime_UTR_variant                |   2428|
|synonymous_variant                 |   1225|
|missense_variant                   |    714|
|5_prime_UTR_variant                |    489|
|splice_region_variant              |    295|
|no_sequence_alteration             |     10|
|stop_gained                        |      6|
|initiator_codon_variant            |      5|
|complex_transcript_variant         |      4|
|splice_donor_variant               |      4|
|coding_sequence_variant            |      2|
|frameshift_variant                 |      2|
|splice_acceptor_variant            |      2|
|stop_lost                          |      2|
|inframe_insertion                  |      1|
|stop_retained_variant              |      1|
:::

::: {.cell-output-display}
Table: Variants by annotation type

|Gene        |  n|
|:-----------|--:|
|Ifi207      | 82|
|Ifi206      | 48|
|Ifi203      | 42|
|Cd244a      | 41|
|Olfr418     | 38|
|Mndal       | 26|
|Ly9         | 23|
|Fcgr2b      | 18|
|Ifi213      | 17|
|Mptx1       | 16|
|Fcrla       | 15|
|Fmn2        | 14|
|Olfr427     | 14|
|Igsf8       | 13|
|Ackr1       | 12|
|Apoa2       | 12|
|Igsf9       | 11|
|Dusp12      | 10|
|Ifi202b     | 10|
|Ifi208      | 10|
|Slamf7      | 10|
|Exo1        |  9|
|Fcgr3       |  8|
|Cd84        |  7|
|Cep170      |  7|
|Itln1       |  7|
|Arhgap30    |  6|
|Cd48        |  6|
|Fcrl6       |  6|
|Ifi205      |  6|
|Olfr1406    |  6|
|Olfr424     |  6|
|Olfr429     |  6|
|Sdccag8     |  6|
|Sdhc        |  6|
|Cfap126     |  5|
|Ifi209      |  5|
|Ifi211      |  5|
|Ccdc190     |  4|
|Copa        |  4|
|Fcer1a      |  4|
|Gm10521     |  4|
|Mpz         |  4|
|Olfr220     |  4|
|Olfr420     |  4|
|Slamf9      |  4|
|Spta1       |  4|
|Aim2        |  3|
|Atf6        |  3|
|Atp1a4      |  3|
|Chml        |  3|
|Gm26620     |  3|
|Gm7694      |  3|
|Hsd17b7     |  3|
|Klhdc9      |  3|
|Ncstn       |  3|
|Nit1        |  3|
|Olfml2b     |  3|
|Olfr218     |  3|
|Olfr248     |  3|
|Olfr417     |  3|
|Olfr432     |  3|
|Slamf8      |  3|
|Uap1        |  3|
|Usp21       |  3|
|Wdr64       |  3|
|Adamts4     |  2|
|Catspere1   |  2|
|Catspere2   |  2|
|Crp         |  2|
|Fcer1g      |  2|
|Fcrlb       |  2|
|Nos1ap      |  2|
|Nuf2        |  2|
|Olfr414     |  2|
|Olfr419     |  2|
|Olfr433     |  2|
|Pfdn2       |  2|
|Rbm8a2      |  2|
|Sh2d1b1     |  2|
|Uhmk1       |  2|
|Alyref2     |  1|
|Fcgr4       |  1|
|Ifi204      |  1|
|Kcnj10      |  1|
|Olfr421-ps1 |  1|
|Pex19       |  1|
|Pigm        |  1|
|Sh2d1b2     |  1|
:::
:::


# Known Associations


Downloaded known total cholesterol variants from https://t2d.hugeamp.org/phenotype.html?phenotype=CHOL


::: {.cell}
::: {.cell-output-display}
Table: High impact variants with human cholesterol associating orthologs

|phenotypes |gene | chromosome| start| end| minP| CHOL:pValue| CHOL:zStat|CHOL:nParam | CHOL:subjects| chiSquared| DB.Class.Key|mouse |
|:----------|:----|----------:|-----:|---:|----:|-----------:|----------:|:-----------|-------------:|----------:|------------:|:-----|
:::

::: {.cell-output-display}
Table: All variants with human cholesterol associating orthologs

|phenotypes |gene | chromosome| start| end| minP| CHOL:pValue| CHOL:zStat|CHOL:nParam | CHOL:subjects| chiSquared| DB.Class.Key|mouse |
|:----------|:----|----------:|-----:|---:|----:|-----------:|----------:|:-----------|-------------:|----------:|------------:|:-----|
:::
:::


# Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}
```
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur ... 10.16

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] stringr_1.5.0 broom_1.0.5   ggplot2_3.4.4 tidyr_1.3.0   dplyr_1.1.3  
[6] readr_2.1.4   knitr_1.44   

loaded via a namespace (and not attached):
 [1] pillar_1.9.0      compiler_4.2.2    tools_4.2.2       bit_4.0.5        
 [5] digest_0.6.33     jsonlite_1.8.7    evaluate_0.22     lifecycle_1.0.3  
 [9] tibble_3.2.1      gtable_0.3.4      pkgconfig_2.0.3   rlang_1.1.1      
[13] cli_3.6.1         rstudioapi_0.15.0 parallel_4.2.2    yaml_2.3.7       
[17] xfun_0.40         fastmap_1.1.1     withr_2.5.1       generics_0.1.3   
[21] vctrs_0.6.4       htmlwidgets_1.6.2 hms_1.1.3         bit64_4.0.5      
[25] grid_4.2.2        tidyselect_1.2.0  glue_1.6.2        R6_2.5.1         
[29] fansi_1.0.5       vroom_1.6.4       rmarkdown_2.25    tzdb_0.4.0       
[33] purrr_1.0.2       magrittr_2.0.3    backports_1.4.1   scales_1.2.1     
[37] htmltools_0.5.6.1 colorspace_2.1-0  utf8_1.2.3        stringi_1.7.12   
[41] munsell_0.5.0     crayon_1.5.2     
```
:::
:::