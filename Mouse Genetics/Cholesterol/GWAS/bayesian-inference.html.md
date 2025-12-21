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
ncd.clump.data <- read_csv(ncd.clump.filename) |>
  filter(gene_biotype == 'protein_coding') |>
  filter(!is.na(external_gene_name))

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


hfd.clump.data <- read_csv(hfd.clump.filename) |>
  filter(gene_biotype == 'protein_coding') |>
  filter(!is.na(external_gene_name))

hfd.clump.summary <-
  hfd.clump.data |>
  group_by(SNP) |>
  summarize(BP = first(BP),
            Min.QTL = first(Min_Position),
            Max.QTL = first(Max_Position),
            Genes = paste(external_gene_name, collapse = ';'),
            Gene.Count = length(ensembl_gene_id)) |>
  rename(Clump=SNP) |>
  mutate(Prior = 1 / Gene.Count)

twas.all.filename <- 'Cholesterol - Expression Associations - Linear Model.csv'
twas.all.data <- read_csv(twas.all.filename)

human.cholesterol.filename <- 'gene_associations_total_cholesterol.csv'
human.chol.data <- read_csv(human.cholesterol.filename)

mgi_file <- 'HOM_MouseHumanSequence.rpt'
#download.file(url=mgi_table,destfile=mgi_file)
mouse_human_genes = read_tsv(mgi_file)

library(tidyr)
mouse.human.table <- mouse_human_genes %>%
  select(`DB Class Key`, `Common Organism Name`,`Symbol`) %>%
  group_by(`DB Class Key`) %>%
  pivot_wider(names_from=`Common Organism Name`,
              values_from=Symbol,
              values_fn=first) %>%
  rename(mouse=`mouse, laboratory`)

human.chol.data.annot <-
  left_join(human.chol.data,
            mouse.human.table,
            by=c("gene"="human")) |>
  arrange(desc(huge)) |>
  distinct(mouse, .keep_all = T)

twas.all.data <- left_join(twas.all.data,
                              human.chol.data.annot |>
                                select(mouse, huge),
                              by=c("symbol"="mouse")) 
```
:::


There is an average of 62.5 genes [12 -155] per NCD QTL and 52.3 genes per HFD QTL [13 -130].  This is filtered to only include protein coding genes with known external gene names.

## Formalizing of observed data

### Coloc data


::: {.cell}

```{.r .cell-code}
coloc.summary.filename <- 'Coloc/Coloc results summary.csv'
coloc.data <- read_csv(coloc.summary.filename)

alpha.coloc <- 1.25

coloc.summary <-
  coloc.data |>
  group_by(Gene, Diet) %>%
  arrange(desc(PP.H4.abf), desc(PP.H3.abf), desc(PP.H2.abf)) %>%
  slice(1) %>%
  ungroup() |>
  mutate(logLR_coloc = alpha.coloc * log(PP.H4.abf/(PP.H2.abf + PP.H3.abf + alpha.coloc))) 
```
:::


From coloc analyses, found in  we have several probabilities:

* $PP.H2$ is the probability of eQTL signal only
* $PP.H3$ is the probability of both traits being associated, but with *different causal* variants
* $PP.H4$ is the probability of both traits being associated and share the *same causal* variant

$LR_{coloc,g}=\frac{P(data∣C_g=0)}{P(data∣C_g=1)}=\alpha \times log(\frac{PP.H4}{PP.H3+PP.H2+\epsilon})$

From this analysis:

- $\epsilon$ is a numerical and conceptual stabilizer, not a biological parameter.  It prevents division by zero, extreme log-odds when posteriors are near zero, and overconfidence from poorly powered regions.  Generally something like $1 \times 10^{-6}$ or more conservatively $1\times 10^{-4}$.
- $\alpha$ is a weighting parameter controlling how strongly coloc evidence influences the final posterior.  Something between $1-1.5$ is reasonable here.  I used 1.25.

### Distance to QTL


::: {.cell}

```{.r .cell-code}
lambda <- 0.1  # relatively weak lambda

ncd.gene.distances <-
  ncd.clump.data %>%
  filter(gene_biotype == "protein_coding") %>%
  mutate(
    Gene_Midpoint = (start_position + end_position) / 2,

    # distance to clump interval (in Mb)
    Distance_to_QTL = case_when(
      Gene_Midpoint >= Min_Position & Gene_Midpoint <= Max_Position ~ 0,
      Gene_Midpoint < Min_Position ~ (Min_Position - Gene_Midpoint) / 1e6,
      Gene_Midpoint > Max_Position ~ (Gene_Midpoint - Max_Position) / 1e6
    ),

    # log-likelihood ratio from distance
    logLR_dist = -lambda * log1p(Distance_to_QTL)
  ) %>%
  select(
    SNP,
    ensembl_gene_id,
    external_gene_name,
    Distance_to_QTL,
    logLR_dist
  )

hfd.gene.distances <-
  hfd.clump.data %>%
  filter(gene_biotype == "protein_coding") %>%
  mutate(
    Gene_Midpoint = (start_position + end_position) / 2,

    # distance to clump interval (in Mb)
    Distance_to_QTL = case_when(
      Gene_Midpoint >= Min_Position & Gene_Midpoint <= Max_Position ~ 0,
      Gene_Midpoint < Min_Position ~ (Min_Position - Gene_Midpoint) / 1e6,
      Gene_Midpoint > Max_Position ~ (Gene_Midpoint - Max_Position) / 1e6
    ),

    # log-likelihood ratio from distance
    logLR_dist = -lambda * log1p(Distance_to_QTL)
  ) %>%
  select(
    SNP,
    ensembl_gene_id,
    external_gene_name,
    Distance_to_QTL,
    logLR_dist
  )
```
:::


Using gene coordinates, we defined a distance-based prior modifier rather than a hard filter. This was implemented as $logLR_{dist,g} = - \lambda \times log(1+d_g)$ where $d_g$ is the distance (in Mb) from gene $g$ to the nearest boundary of the LD block containing the QTL, and $\lambda$ controls the influence of spatial proximity on the final posterior.  $\lambda$ was assigned a value of 0.1, reflecting a deliberately weak spatial prior.

Genes overlapping the QTL interval were assigned $\log LR_{distg,d}=0$, reflecting the limited positional resolution of DO mouse QTLs. The distance penalty was intentionally weak, such that genes located several megabases outside the interval incurred only modest log-odds penalties, allowing strong genetic, transcriptional, or colocalization evidence to outweigh spatial proximity.  This formulation reflects the broad haplotype structure of DO mouse populations, in which association signals often span multiple megabases and peak position provides limited additional localization within the QTL interval.

### Gene Expression Proxy

Gene expression–trait associations were treated as supportive but non-causal evidence, conditional on genetic localization. For each gene ($g$), we defined a standardized association statistic $Z_{expr, g}=\frac{\beta_g}{SE_g}$ where $\beta_g$ is the beta coefficient for the sex-adjusted scaled cholesterol-gene expression linear model and $SE_g$.

We then defined an expression-based log-likelihood ratio component as $\log LR_{expr,g}=\omega \times \frac{Z_{expr,g}^2}{2}$ which corresponds to the log-likelihood ratio for testing $\beta_g=0$ under the normal a normal approximation. The contribution was capped at a maximum value to prevent gene expression evidence from dominating the posterior. This formulation rewards strong expression–trait associations while avoiding penalization of weak or noisy effects and reflects the non-causal but informative role of expression evidence in the model.


::: {.cell}

```{.r .cell-code}
omega <- 0.2
max_reward <- 1.0

moduleB_logLR <- function(beta, se,
                          omega = 0.2,
                          max_reward = 1.0) {

  Z <- beta / se
  Z[!is.finite(Z)] <- 0  # handles NA, Inf, -Inf

  logLR <- omega * (Z^2) / 2
  logLR <- pmin(logLR, max_reward)

  return(logLR)
}
```
:::


The weighting parameters included $\omega$=0.2, reflecting the deliberately low weight assigned to expression–trait associations, and a maximum reward of 1, which prevents this non-causal evidence from dominating the posterior. As a result, expression associations can support—but not override—genetic evidence when nominating causal genes.

## Human Data Proxy

Bayes factors for associations with total cholesterol are downloaded from the [Common Metabolic Disease Knowledge Portal](https://hugeamp.org/).  This is parameterized as: 

$$logLR_{human,g} = \{ \frac{log(BF_{human,g})}{0} \frac{\text{if it exists}}{\text{if it does not exist}}$$

# Analyses

## Chromosome 5 QTL (HFD specific)

This was first for the '5_117508066_B_E' clump, which is centered on the 123629774 bp position.


::: {.cell}

```{.r .cell-code}
twas.all.data %>%
  right_join(hfd.clump.data %>%
               distinct(ensembl_gene_id, .keep_all = T) |>
               filter(SNP=="5_117508066_B_E") |> #clump boundaries for the 123629774 centered QTL
               select(SNP,ensembl_gene_id, external_gene_name, CHR, start_position, end_position), by=c("ID"="ensembl_gene_id")) |> 
  mutate(logLR_expression = moduleB_logLR(beta, se)) |>
  select(SNP,chr,symbol,logLR_expression,beta,huge) |>
  filter(!is.na(symbol)) |>
  left_join(hfd.clump.summary |>
              select(Clump,Prior,Gene.Count),by=c("SNP"="Clump")) |> 
  inner_join(coloc.summary |> select(Gene,Diet, logLR_coloc), by=c("symbol"="Gene")) |> 
  left_join(hfd.gene.distances |> 
              select(logLR_dist,external_gene_name),
            by=c("symbol"="external_gene_name")) |>
  filter(!is.na(SNP)) |>
  mutate(logLR_prior = log(Prior)) |>
  mutate(logLR_human = ifelse(!is.na(huge), log(huge), 0)) |>
  mutate(log_odds = logLR_prior + logLR_coloc + logLR_dist + logLR_expression + logLR_human) |>
    mutate(posterior = plogis(log_odds)) |>
  select(SNP,symbol, posterior, log_odds, Prior, logLR_prior, logLR_coloc, logLR_dist, beta, logLR_expression, logLR_human) |>
  arrange(-log_odds) ->
  chr5.summary.data

kable(chr5.summary.data, caption="liklihood of each gene in the chromosome 5 QTL")
```

::: {.cell-output-display}


Table: liklihood of each gene in the chromosome 5 QTL

|SNP             |symbol | posterior|  log_odds|     Prior| logLR_prior| logLR_coloc| logLR_dist|      beta| logLR_expression| logLR_human|
|:---------------|:------|---------:|---------:|---------:|-----------:|-----------:|----------:|---------:|----------------:|-----------:|
|5_117508066_B_E |Mvk    | 0.0180006| -3.999186| 0.0128205|   -4.356709|  -0.6424775|          0| -16.14079|                1|           0|
|5_117508066_B_E |Mvk    | 0.0169676| -4.059338| 0.0128205|   -4.356709|  -0.7026293|          0| -16.14079|                1|           0|


:::
:::


This was first for the '5_117508066_B_E' clump, which is centered on the 123629774 bp position.


::: {.cell}

```{.r .cell-code}
twas.all.data %>%
  right_join(hfd.clump.data %>%
               distinct(ensembl_gene_id, .keep_all = T) |>
               filter(SNP=="5_123629774_B_E") |> #clump boundaries for the 123629774 centered QTL
               select(SNP,ensembl_gene_id, external_gene_name, CHR, start_position, end_position), by=c("ID"="ensembl_gene_id")) |> 
  mutate(logLR_expression = moduleB_logLR(beta, se)) |>
  select(SNP,chr,symbol,logLR_expression,beta,huge) |>
  filter(!is.na(symbol)) |>
  left_join(hfd.clump.summary |>
              select(Clump,Prior,Gene.Count),by=c("SNP"="Clump")) |> 
  inner_join(coloc.summary |> select(Gene,Diet, logLR_coloc), by=c("symbol"="Gene")) |> 
  left_join(hfd.gene.distances |> 
              select(logLR_dist,external_gene_name),
            by=c("symbol"="external_gene_name")) |>
  filter(!is.na(SNP)) |>
  mutate(logLR_prior = log(Prior)) |>
  mutate(logLR_human = ifelse(!is.na(huge), log(huge), 0)) |>
  mutate(log_odds = logLR_prior + logLR_coloc + logLR_dist + logLR_expression + logLR_human) |>
    mutate(posterior = plogis(log_odds)) |>
  select(SNP,symbol, posterior, log_odds, Prior, logLR_prior, logLR_coloc, logLR_dist, beta, logLR_expression, logLR_human) |>
  arrange(-log_odds) ->
  chr5.summary.data

kable(chr5.summary.data, caption="liklihood of each gene in the second chromosome 5 QTL")
```

::: {.cell-output-display}


Table: liklihood of each gene in the second chromosome 5 QTL

|SNP             |symbol | posterior|  log_odds|     Prior| logLR_prior| logLR_coloc| logLR_dist|      beta| logLR_expression| logLR_human|
|:---------------|:------|---------:|---------:|---------:|-----------:|-----------:|----------:|---------:|----------------:|-----------:|
|5_123629774_B_E |Scarb1 | 0.2062124| -1.347909| 0.0083333|   -4.787492|  -0.3944073|          0| -1.314871|        0.0157955|    3.818194|
|5_123629774_B_E |Scarb1 | 0.2062124| -1.347909| 0.0083333|   -4.787492|  -0.3944073|          0| -1.314871|        0.0157955|    3.818194|


:::
:::

## Chromosome 13 QTL (HFD specific)


::: {.cell}

```{.r .cell-code}
twas.all.data %>%
  right_join(hfd.clump.data %>%
               filter(CHR==13) %>%
               select(SNP,ensembl_gene_id, external_gene_name, CHR, start_position, end_position),
             by = c('ID'='ensembl_gene_id')) %>%
  mutate(logLR_expression = moduleB_logLR(beta, se)) |>
  select(SNP,chr,symbol,logLR_expression,beta,huge) |>
  filter(!is.na(symbol)) |>
  left_join(hfd.clump.summary |>
              select(Clump,Prior,Gene.Count),by=c("SNP"="Clump")) |>
  full_join(coloc.summary |> select(Gene,Diet, logLR_coloc), by=c("symbol"="Gene")) |>
  left_join(hfd.gene.distances |> select(logLR_dist,external_gene_name), by=c("symbol"="external_gene_name")) |>
  filter(!is.na(SNP)) |>
    mutate(logLR_prior = log(Prior)) |>
  mutate(logLR_human = ifelse(!is.na(huge), log(huge), 0)) |>
  mutate(log_odds = logLR_prior + logLR_coloc + logLR_dist + logLR_expression + logLR_human) |>
    mutate(posterior = plogis(log_odds)) |>
  select(SNP,symbol, posterior, log_odds, Prior, logLR_prior, logLR_coloc, logLR_dist, beta, logLR_expression, logLR_human) |>
  arrange(-logLR_expression)  ->
  chr13.summary.data

kable(chr13.summary.data, caption="liklihood of each gene in the chromosome 13 QTL")
```

::: {.cell-output-display}


Table: liklihood of each gene in the chromosome 13 QTL

|SNP             |symbol    | posterior|  log_odds|     Prior| logLR_prior| logLR_coloc| logLR_dist|        beta| logLR_expression| logLR_human|
|:---------------|:---------|---------:|---------:|---------:|-----------:|-----------:|----------:|-----------:|----------------:|-----------:|
|13_30180778_E_C |Cdkal1    |        NA|        NA| 0.0153846|   -4.174387|          NA|          0| -11.0765308|        1.0000000|    0.000000|
|13_30180778_E_C |Wrnip1    |        NA|        NA| 0.0153846|   -4.174387|          NA|          0| -18.5230649|        1.0000000|    0.000000|
|13_30180778_E_C |Foxq1     |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -5.2134928|        1.0000000|    0.000000|
|13_30180778_E_C |Serpinb6b |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|   9.8059085|        1.0000000|    0.000000|
|13_30180778_E_C |Serpinb1a |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -4.1688716|        1.0000000|    0.000000|
|13_30180778_E_C |Agtr1a    |        NA|        NA| 0.0153846|   -4.174387|          NA|          0| -15.8968520|        1.0000000|    0.000000|
|13_30180778_E_C |Dcdc2a    |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -5.4837801|        0.9326417|    0.000000|
|13_30180778_E_C |Ripk1     |        NA|        NA| 0.0153846|   -4.174387|          NA|          0| -10.0734756|        0.7815369|    0.000000|
|13_30180778_E_C |Bphl      |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -7.4225183|        0.5469114|    0.000000|
|13_30180778_E_C |Uqcrfs1   |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -8.3692532|        0.5070385|    0.000000|
|13_30180778_E_C |Nqo2      |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -6.1245074|        0.4639308|    0.000000|
|13_30180778_E_C |Exoc2     |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -7.6802585|        0.3203858|    0.000000|
|13_30180778_E_C |Dusp22    |         0| -23.50833| 0.0153846|   -4.174387|   -19.65112|          0|   6.4069483|        0.3171717|    0.000000|
|13_30180778_E_C |E2f3      |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|   5.2752424|        0.2544446|    3.806662|
|13_30180778_E_C |Serpinb6a |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -4.5508752|        0.2254086|    0.000000|
|13_30180778_E_C |Slc22a23  |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -4.9213482|        0.1553026|    0.000000|
|13_30180778_E_C |Serpinb9  |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -3.2339024|        0.1313383|    0.000000|
|13_30180778_E_C |Tubb2a    |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -0.9144915|        0.0647426|    0.000000|
|13_30180778_E_C |Prl8a1    |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -1.6307464|        0.0532674|    0.000000|
|13_30180778_E_C |Psmg4     |        NA|        NA| 0.0153846|   -4.174387|          NA|          0|  -1.4098899|        0.0294358|    0.000000|


:::
:::


## Chromosome 18 QTL (NCD specific)


::: {.cell}

```{.r .cell-code}
twas.all.data %>%
  right_join(ncd.clump.data %>%
               filter(CHR==18) %>%
               select(SNP,ensembl_gene_id, external_gene_name, CHR, start_position, end_position),
             by = c('ID'='ensembl_gene_id')) %>%
  mutate(logLR_expression = moduleB_logLR(beta, se)) |>
  select(SNP,chr,symbol,logLR_expression,beta,huge) |>
  filter(!is.na(symbol)) |>
  left_join(ncd.clump.summary |>
              select(Clump,Prior,Gene.Count),by=c("SNP"="Clump")) |>
  full_join(coloc.summary |> select(Gene,Diet, logLR_coloc), by=c("symbol"="Gene")) |>
  left_join(ncd.gene.distances |> select(logLR_dist,external_gene_name), by=c("symbol"="external_gene_name")) |>
  filter(!is.na(SNP)) |>
    mutate(logLR_prior = log(Prior)) |>
  mutate(logLR_human = ifelse(!is.na(huge), log(huge), 0)) |>
  mutate(log_odds = logLR_prior + logLR_coloc + logLR_dist + logLR_expression + logLR_human) |>
    mutate(posterior = plogis(log_odds)) |>
  select(SNP,symbol, posterior, log_odds, Prior, logLR_prior, logLR_coloc, logLR_dist, beta, logLR_expression, logLR_human) |>
  arrange(-log_odds)  ->
  chr18.summary.data

kable(chr18.summary.data, caption="liklihood of each gene in the chromosome 18 QTL")
```

::: {.cell-output-display}


Table: liklihood of each gene in the chromosome 18 QTL

|SNP             |symbol  | posterior|   log_odds|     Prior| logLR_prior| logLR_coloc| logLR_dist|       beta| logLR_expression| logLR_human|
|:---------------|:-------|---------:|----------:|---------:|-----------:|-----------:|----------:|----------:|----------------:|-----------:|
|18_46410922_G_A |Dmxl1   | 0.0334071|  -3.365010| 0.0227273|    -3.78419|  -0.5808205|   0.000000| -16.540550|        1.0000000|           0|
|18_46410922_G_A |Cdo1    | 0.0326668|  -3.388185| 0.0227273|    -3.78419|  -0.6039950|   0.000000| -10.423712|        1.0000000|           0|
|18_46410922_G_A |Mcc     | 0.0243427|  -3.690879| 0.0227273|    -3.78419|  -0.9066892|   0.000000| -14.665844|        1.0000000|           0|
|18_46410922_G_A |Ythdc2  | 0.0083569|  -4.776278| 0.0227273|    -3.78419|  -1.5386081|   0.000000|  -7.773632|        0.5465201|           0|
|18_46410922_G_A |Kcnn2   | 0.0028732|  -5.849440| 0.0227273|    -3.78419|  -3.0652501|   0.000000| -17.365697|        1.0000000|           0|
|18_46410922_G_A |Ap3s1   | 0.0001126|  -9.091765| 0.0227273|    -3.78419|  -5.5987425|   0.000000|   5.295224|        0.2911672|           0|
|18_46410922_G_A |Dcp2    | 0.0000576|  -9.762419| 0.0227273|    -3.78419|  -6.5178999|   0.000000|  -9.187939|        0.5396704|           0|
|18_46410922_G_A |Tcerg1  | 0.0000436| -10.039880| 0.0227273|    -3.78419|  -7.2556904|   0.000000| -12.041842|        1.0000000|           0|
|18_46410922_G_A |Tcerg1  | 0.0000397| -10.133446| 0.0227273|    -3.78419|  -7.3492565|   0.000000| -12.041842|        1.0000000|           0|
|18_46410922_G_A |Fem1c   | 0.0000207| -10.786933| 0.0227273|    -3.78419|  -7.7043111|   0.000000| -10.311275|        0.7015675|           0|
|18_46410922_G_A |Sema6a  |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -4.388689|        0.2057660|           0|
|18_46410922_G_A |Pggt1b  |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -8.989326|        0.7571650|           0|
|18_46410922_G_A |Dtwd2   |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -2.762396|        0.0662733|           0|
|18_46410922_G_A |Atg12   |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -5.657356|        0.2484670|           0|
|18_46410922_G_A |Tmed7   |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -5.326377|        0.2239222|           0|
|18_46410922_G_A |Commd10 |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -4.699061|        0.2983901|           0|
|18_46410922_G_A |Eif3j2  |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -4.670281|        0.2128359|           0|
|18_46410922_G_A |Eif1a   |        NA|         NA| 0.0227273|    -3.78419|          NA|   0.000000|  -2.435114|        0.0960097|           0|
|18_46410922_G_A |Tnfaip8 |        NA|         NA| 0.0227273|    -3.78419|          NA|  -0.000667|  -2.217297|        0.0958360|           0|


:::
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
Running under: macOS Tahoe 26.2

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


