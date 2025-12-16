---
title: "Bayesian inference to identify locus-specific causal genes for cholesterol QTLs"
author: "Dave Bridges"
date: "2025-12-15"
editor: source
format: 
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
    fig-path: "figures/"
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


# Purpose

To set up and execute a causal framework to evaluate causal genes in individual QTLs we have identified for cholesterol.

# Framework

This used chatGPT to help draft this code.  The presumption is that we have several pieces of information

* The presence of a cholesterol QTL either diet specific or not
* Colocalization analyses using coloc/SuSiE

Our goal is to compute $P(Cg​=1∣\text{all observed data})$

Our *prior probablity* for each gene sparse, and is based on $P(Cg​=1)=\pi_o\approx\frac{1}{N_\text{genes in locus}}$


::: {.cell}

```{.r .cell-code}
ncd.clump.filename <- 'ld-calculations/Annotated Gene Clumps NCD QTLs.csv'
hfd.clump.filename <- 'ld-calculations/Annotated Gene Clumps HFD QTLs.csv'

library(readr)
ncd.clump.data <- read_csv(ncd.clump.filename) 

ncd.clump.summary <-
  ncd.clump.data |>
  group_by(SNP) |>
  filter(gene_biotype == 'protein_coding') |>
  summarize(BP = first(BP),
            Min.QTL = first(Min_Position),
            Max.QTL = first(Max_Position),
            Genes = paste(external_gene_name, collapse = ';'),
            Gene.Count = length(ensembl_gene_id)) |>
  rename(Clump=SNP) |>
  mutate(Prior = 1 / Gene.Count)


hfd.clump.data <- read_csv(hfd.clump.filename) 

hfd.clump.summary <-
  hfd.clump.data |>
  group_by(SNP) |>
  filter(gene_biotype == 'protein_coding') |>
  summarize(BP = first(BP),
            Min.QTL = first(Min_Position),
            Max.QTL = first(Max_Position),
            Genes = paste(external_gene_name, collapse = ';'),
            Gene.Count = length(ensembl_gene_id)) |>
  rename(Clump=SNP) |>
  mutate(Prior = 1 / Gene.Count)

twas.all.filename <- 'Cholesterol - Expression Associations - Linear Model.csv'
twas.all.data <- read_csv(twas.all.filename)
```
:::


There is an average of 62.7 genes [12 -156] per NCD QTL and 52.5 genes per HFD QTL [13 -130] .

## Formalizing of observed data

### Coloc data


::: {.cell}

```{.r .cell-code}
coloc.summary.filename <- 'Coloc/Coloc results summary.csv'
coloc.data <- read_csv(coloc.summary.filename)

alpha.coloc <- 1

coloc.summary <-
  coloc.data |>
  group_by(Gene, Diet) %>%
  arrange(desc(PP.H4.abf), desc(PP.H3.abf), desc(PP.H2.abf)) %>%
  slice(1) %>%
  ungroup() |>
  mutate(logLR.coloc = alpha.coloc * log(PP.H4.abf/(PP.H2.abf + PP.H3.abf + alpha.coloc))) 
```
:::


From coloc analyses, found in  we have several probabilities:

* $PP.H2$ is the probability of eQTL signal only
* $PP.H3$ is the probability of both traits being associated, but with *different causal* variants
* $PP.H4$ is the probability of both traits being associated and share the *same causal* variant

$LR_{coloc,g}=\frac{P(data∣C_g=0)}{P(data∣C_g=1)}=\alpha \times log(\frac{PP.H4}{PP.H3+PP.H2+\epsilon})$

From this analysis:

- $\epsilon$ is a numerical and conceptual stabilizer, not a biological parameter.  It prevents division by zero, extreme log-odds when posteriors are near zero, and overconfidence from poorly powered regions.  Generally something like $1 \times 10^{-6}$ or more conservatively $1\times 10^{-4}$.
- $\alpha$ is a weighting parameter controlling how strongly coloc evidence influences the final posterior.  Something between $1-1.5$ is reasonable here.  I used 1.

### Distance to QTL

Using gene coordinates to define a distance-based prior modifier, rather than a filter.  This will be defined as $log(LR_{dist,g}λ⋅lo\lambda1+dg​)$_ where $\l4ambda$ is a scaling parameter controlling the influence of distance on the final posterior.  $d_g$ is the distance from the gene to the lead SNP for the QTL.


::: {.cell}

```{.r .cell-code}
lambda <- 1.0
ncd.clump.data %>%
  mutate(Gene_Midpoint = (start_position + end_position) / 2) |>
  mutate(Distance_to_QTL = abs(Gene_Midpoint - BP)/1E6) |> #in Mb
  mutate(logLR_dist = -lambda * log(1+Distance_to_QTL)) |>
  select(SNP,ensembl_gene_id, external_gene_name, Distance_to_QTL, logLR_dist) ->
  ncd.gene.distances
```
:::


We incorporated a weak distance-based prior that modestly favored genes proximal to the QTL peak while allowing strong genetic or transcriptional evidence to outweigh spatial proximity. The distance penalty was calibrated such that a gene 1 Mb from the peak incurred less than `one1 log-odds unit of penalty.

## Gene Expression Proxy

This is not causal by itself, but informative conditional on genetics.

To do this, first we defined a standardized effect: $Z_{expr,g} = \frac{β_{g}}{SE_{g}}$ where $β_{g}$ is the effect estimate from the linear model of expression vs cholesterol and $SE_{g}$ is the standard error of that estimate.


::: {.cell}

```{.r .cell-code}
mu_med    <- 1.0
sigma_med <- 1.0

moduleB_logLR <- function(beta, se,
                          max_reward = 1.0,
                          weight = 0.2) {

  Z <- beta / se

  # reward strong associations, but do not penalize weak ones
  logLR <- weight * (Z^2) / 2

  # cap the contribution
  logLR <- pmin(logLR, max_reward)

  return(logLR)
}

twas.all.data %>%
  right_join(hfd.clump.data %>%
               select(SNP,ensembl_gene_id, external_gene_name, CHR, start_position, end_position) %>%
               distinct(),
             by = c('ID'='ensembl_gene_id')) %>%
  mutate(logLR.expression = moduleB_logLR(beta, se)) |>
  select(chr,symbol,logLR.expression) ->
  hfd.gene.expression.proxy
```
:::


Then define a likelihood ratio component as $LR_{expr,g} = logit^{-1}(Z_{expr,g})$.

# Analyses

## Chromosome 18 QTL NCD specific


::: {.cell}

```{.r .cell-code}
twas.all.data %>%
  right_join(ncd.clump.data %>%
               filter(CHR==18) %>%
               select(SNP,ensembl_gene_id, external_gene_name, CHR, start_position, end_position),
             by = c('ID'='ensembl_gene_id')) %>%
  mutate(logLR_expression = moduleB_logLR(beta, se)) |>
  select(SNP,chr,symbol,logLR_expression) |>
  filter(!is.na(symbol)) |>
  left_join(ncd.clump.summary |>
              select(Clump,Prior,Gene.Count),by=c("SNP"="Clump")) |>
  full_join(coloc.summary |> select(Gene,Diet, logLR.coloc), by=c("symbol"="Gene")) |>
  left_join(ncd.gene.distances |> select(logLR_dist,external_gene_name), by=c("symbol"="external_gene_name")) |>
  filter(!is.na(SNP)) |>
  mutate(log_prior = log(Prior)) |>
  mutate(log_odds = log_prior + logLR.coloc + logLR_dist + logLR_expression) |>
  select(symbol, log_odds, log_prior, logLR.coloc, logLR_dist, logLR_expression) |>
  mutate(posterior = plogis(log_odds)) |>
  arrange(-log_odds) ->
  chr18.summary.data

kable(chr18.summary.data, caption="liklihood of each gene in the chromosome 18 QTL")
```

::: {.cell-output-display}


Table: liklihood of each gene in the chromosome 18 QTL

|symbol  |   log_odds| log_prior| logLR.coloc| logLR_dist| logLR_expression| posterior|
|:-------|----------:|---------:|-----------:|----------:|----------------:|---------:|
|Cdo1    |  -3.436049|  -3.78419|  -0.2851774| -0.3666820|        1.0000000| 0.0311876|
|Mcc     |  -4.307368|  -3.78419|  -0.5451403| -0.9780381|        1.0000000| 0.0132900|
|Dmxl1   |  -4.579920|  -3.78419|  -0.2650722| -1.5306580|        1.0000000| 0.0101516|
|Ythdc2  |  -5.196990|  -3.78419|  -1.0759163| -0.8834038|        0.5465201| 0.0055027|
|Kcnn2   |  -5.695841|  -3.78419|  -2.3235957| -0.5880552|        1.0000000| 0.0033487|
|Ap3s1   |  -8.250580|  -3.78419|  -4.3597919| -0.3977662|        0.2911672| 0.0002610|
|Fem1c   |  -9.326097|  -3.78419|  -6.0315915| -0.2118832|        0.7015675| 0.0000891|
|Dcp2    |  -9.362995|  -3.78419|  -5.0623811| -1.0560942|        0.5396704| 0.0000858|
|Tcerg1  | -10.014247|  -3.78419|  -5.6752411| -1.5548163|        1.0000000| 0.0000448|
|Tcerg1  | -10.089125|  -3.78419|  -5.7501185| -1.5548163|        1.0000000| 0.0000415|
|Sema6a  |         NA|  -3.78419|          NA| -0.7052630|        0.2057660|        NA|
|Pggt1b  |         NA|  -3.78419|          NA| -0.0194918|        0.7571650|        NA|
|Dtwd2   |         NA|  -3.78419|          NA| -1.4924586|        0.0662733|        NA|
|Atg12   |         NA|  -3.78419|          NA| -0.3770704|        0.2484670|        NA|
|Tmed7   |         NA|  -3.78419|          NA| -0.2631563|        0.2239222|        NA|
|Commd10 |         NA|  -3.78419|          NA| -0.5981192|        0.2983901|        NA|
|Eif3j2  |         NA|  -3.78419|          NA| -1.3353300|        0.2128359|        NA|
|Eif1a   |         NA|  -3.78419|          NA| -0.2845949|        0.0960097|        NA|
|Tnfaip8 |         NA|  -3.78419|          NA| -1.5613909|        0.0958360|        NA|
|Gm3650  |         NA|  -3.78419|          NA| -1.3360175|        0.2165430|        NA|


:::

```{.r .cell-code}
#df$log_odds <-
#log_prior +
#df$logLR_coloc +
#df$logLR_B +
#df$logLR_dist +
#df$logLR_human
```
:::


# Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Detroit
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] knitr_1.50      lubridate_1.9.4 forcats_1.0.1   stringr_1.6.0  
 [5] dplyr_1.1.4     purrr_1.2.0     readr_2.1.6     tidyr_1.3.1    
 [9] tibble_3.3.0    ggplot2_4.0.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] bit_4.6.0          gtable_0.3.6       jsonlite_2.0.0     crayon_1.5.3      
 [5] compiler_4.5.2     tidyselect_1.2.1   parallel_4.5.2     dichromat_2.0-0.1 
 [9] scales_1.4.0       yaml_2.3.12        fastmap_1.2.0      R6_2.6.1          
[13] generics_0.1.4     htmlwidgets_1.6.4  pillar_1.11.1      RColorBrewer_1.1-3
[17] tzdb_0.5.0         rlang_1.1.6        stringi_1.8.7      xfun_0.54         
[21] S7_0.2.1           bit64_4.6.0-1      timechange_0.3.0   cli_3.6.5         
[25] withr_3.0.2        magrittr_2.0.4     digest_0.6.39      grid_4.5.2        
[29] vroom_1.6.7        rstudioapi_0.17.1  hms_1.1.4          lifecycle_1.0.4   
[33] vctrs_0.6.5        evaluate_1.0.5     glue_1.8.0         farver_2.1.2      
[37] rmarkdown_2.30     tools_4.5.2        pkgconfig_2.0.3    htmltools_0.5.9   
```


:::
:::


