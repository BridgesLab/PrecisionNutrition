---
title: "Liver mRNA TWAS Analysis"
author: "Dave Bridges"
editor: source
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
---


::: {.cell}

```{.r .cell-code}
library(readr)
library(dplyr)
library(tidyr)
library(broom)
library(tidyverse)

rna.datatfile <- 'dataset.mrna.Svenson_DO_HFD.v12.Rds'
all.data <- readRDS(rna.datatfile)

rna.data.raw <- all.data$data$raw %>% as.data.frame()

library(janitor)
rna.data.raw %>% t %>% as.data.frame -> rna.data

gene.annotation <- all.data$annot.mrna %>%
  mutate(start.true = start*1E6,
         end.true = end*1E6)
```
:::


Pulled out phenotype data


::: {.cell}

```{.r .cell-code}
phenotype.datatfile <- 'dataset.phenotype.Svenson_DO_HFD.v12.Rds'
phenotype.data <- readRDS(phenotype.datatfile)

library(tibble)
chol.phenotypes <- 
  phenotype.data$data$raw %>% 
  as.data.frame %>%
  select(chol2) %>%
  mutate(chol2.scale = scale(chol2,center=T)) %>%
  rownames_to_column("ID")

library(forcats)
group.data <- 
  phenotype.data$annot.samples %>%
  select(mouse.id,sex,diet) %>% 
  column_to_rownames('mouse.id') %>% 
  mutate_all(as.factor) %>%
  rownames_to_column('ID')

pheno.data.all <- full_join(chol.phenotypes,group.data,by="ID")
```
:::


## Plots for Cholesterol-Liver mRNA Expression


::: {.cell}

```{.r .cell-code}
library(tibble)
genes.of.interest <- c('Chrna2','Srebf1','Clstn3','Kif13b','Ascc3','Apoa1','Apoa2','Lpin1','Cdkal1','Bdh1', 'Hmgcs2', 'Hmgcl','Hmgcr','Apoa2','Hnf1a','Mmab','Ldlr','Slc16a1','Slc16a6','Slc16a7','Mvk') #error with Slc16a3/8

for (gene in genes.of.interest) {

  gene.id <- filter(gene.annotation, symbol==gene) %>% pull(gene.id)
gene.data <- 
  rna.data[gene.id,] %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column('ID') %>% 
  left_join(pheno.data.all)

ggplot(gene.data %>% filter(!(is.na(sex))),
       aes(y=chol2,x=get(gene.id))) +
  geom_point(aes(col=sex)) +
  stat_smooth(method="lm") +
  labs(title=gene,
       x="Liver mRNA Expression",
       y="Cholesterol (mg/dL)") +
  theme_classic(base_size=16) -> plot.all
print(plot.all)

ggsave(filename=paste0("figures/cholesterol-",gene,"-plot.pdf"),
       device=c("pdf"))
ggsave(filename=paste0("figures/cholesterol-",gene,"-plot.png"),
       device=c("png"))

ggplot(gene.data %>% filter(!(is.na(sex))),
       aes(y=chol2,x=get(gene.id))) +
  geom_point(aes(col=diet)) +
  facet_grid(~sex) +
  stat_smooth(method="lm") +
  labs(title=gene,
       x="Liver mRNA Expression",
       y="Cholesterol (mg/dL)") +
  theme_classic(base_size=16) ->
  plot.sex
  
print(plot.sex)
ggsave(filename=paste0("figures/cholesterol-sex-",gene,"-plot.pdf"),
       device=c("pdf"))
ggsave(filename=paste0("figures/cholesterol-sex-",gene,"-plot.png"),
       device=c("png"))



ggplot(gene.data %>% filter(!(is.na(sex))),
       aes(y=chol2,x=get(gene.id))) +
  geom_point(aes(col=sex)) +
  facet_grid(~diet) +
  stat_smooth(method="lm") +
  labs(title=gene,
       x="Liver mRNA Expression",
       y="Cholesterol (mg/dL)") +
  theme_classic(base_size=16) ->
  plot.sex
  
print(plot.sex)
ggsave(filename=paste0("figures/cholesterol-diet-",gene,"-plot.pdf"),
       device=c("pdf"))
ggsave(filename=paste0("figures/cholesterol-diet-",gene,"-plot.png"),
       device=c("png"))
}
```

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-1.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-2.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-3.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-4.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-5.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-6.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-7.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-8.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-9.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-10.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-11.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-12.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-13.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-14.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-15.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-16.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-17.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-18.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-19.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-20.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-21.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-22.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-23.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-24.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-25.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-26.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-27.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-28.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-29.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-30.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-31.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-32.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-33.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-34.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-35.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-36.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-37.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-38.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-39.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-40.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-41.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-42.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-43.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-44.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-45.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-46.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-47.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-48.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-49.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-50.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-51.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-52.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-53.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-54.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-55.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-56.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-57.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-58.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-59.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-60.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-61.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-62.png){width=672}
:::

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/plots-63.png){width=672}
:::
:::


## Summary Tables

Loaded in the results files from the linear model analyses (diet and sex specific)


::: {.cell}

```{.r .cell-code}
all.datafile <- 'Cholesterol - Expression Associations - Linear Model.csv'
ncd.datafile <- 'Cholesterol - Expression Associations - NCD Mice - Linear Model.csv'
hfd.datafile <- 'Cholesterol - Expression Associations - HFD Mice - Linear Model.csv'
male.datafile <- 'Cholesterol - Expression Associations - Male Mice - Linear Model.csv'
female.datafile <- 'Cholesterol - Expression Associations - Female Mice - Linear Model.csv'

all.data <- read_csv(all.datafile)
ncd.data <- read_csv(ncd.datafile)
hfd.data <- read_csv(hfd.datafile)
male.data <- read_csv(male.datafile)
female.data <- read_csv(female.datafile)
```
:::


The datafiles were:

* All Data - Cholesterol - Expression Associations - Linear Model.csv, this is adjusted for sex
* NCD Data - Cholesterol - Expression Associations - NCD Mice - Linear Model.csv, this is adjustd for sex
* HFD Data = Cholesterol - Expression Associations - HFD Mice - Linear Model.csv, this is adjusted for sex
* Male Data = Cholesterol - Expression Associations - Male Mice - Linear Model.csv, this is adjusted for diet
* Female Data = Cholesterol - Expression Associations - Female Mice - Linear Model.csv, this is adjusted for diet


::: {.cell}

```{.r .cell-code}
column.order <- c("All","NCD","HFD","Males","Females") 
library(knitr)
bind_rows(all.data %>% filter(symbol %in% genes.of.interest) %>% mutate(model="All"),
          ncd.data %>% filter(symbol %in% genes.of.interest) %>% mutate(model="NCD"),
          hfd.data %>% filter(symbol %in% genes.of.interest) %>% mutate(model="HFD"),
          male.data %>% filter(symbol %in% genes.of.interest) %>% mutate(model="Males"),
          female.data %>% filter(symbol %in% genes.of.interest) %>% mutate(model="Females")) %>%
  arrange(symbol,model) %>%
  mutate(model = fct_relevel(model, column.order)) %>%
  relocate(symbol) -> summary.data

summary.data %>%
  kable(caption="Summary of effect estimates for genes of interest in different models.")
```

::: {.cell-output-display}
Table: Summary of effect estimates for genes of interest in different models.

|symbol  |ID                 |        beta|       se|      pval|         r2|      padj|chr | start.true|  end.true|model   |
|:-------|:------------------|-----------:|--------:|---------:|----------:|---------:|:---|----------:|---------:|:-------|
|Apoa1   |ENSMUSG00000032083 |   5.6695008| 2.588110| 0.0289726|  0.0080176| 0.0676468|9   |   46228580|  46230466|All     |
|Apoa1   |ENSMUSG00000032083 |   3.0167813| 2.832740| 0.2879689|  0.1254920| 0.9998133|9   |   46228580|  46230466|Females |
|Apoa1   |ENSMUSG00000032083 |  -5.1453763| 3.268727| 0.1168519|  0.0642907| 0.1743814|9   |   46228580|  46230466|HFD     |
|Apoa1   |ENSMUSG00000032083 |   3.0167813| 2.832740| 0.2879689|  0.1254920| 0.9998133|9   |   46228580|  46230466|Males   |
|Apoa1   |ENSMUSG00000032083 |   3.0167813| 2.832740| 0.2879689|  0.1254920| 0.9998133|9   |   46228580|  46230466|NCD     |
|Apoa2   |ENSMUSG00000005681 |  -2.1914572| 2.581145| 0.3963004| -0.0005943| 0.5090534|1   |  171225054| 171226379|All     |
|Apoa2   |ENSMUSG00000005681 |   0.7837223| 2.775690| 0.7779188|  0.1216188| 0.9998133|1   |  171225054| 171226379|Females |
|Apoa2   |ENSMUSG00000005681 |  -4.8193875| 3.459749| 0.1649862|  0.0620941| 0.2279211|1   |  171225054| 171226379|HFD     |
|Apoa2   |ENSMUSG00000005681 |   0.7837223| 2.775690| 0.7779188|  0.1216188| 0.9998133|1   |  171225054| 171226379|Males   |
|Apoa2   |ENSMUSG00000005681 |   0.7837223| 2.775690| 0.7779188|  0.1216188| 0.9998133|1   |  171225054| 171226379|NCD     |
|Ascc3   |ENSMUSG00000038774 | -10.4271962| 3.629316| 0.0042497|  0.0152003| 0.0150613|10  |   50592669|  50851202|All     |
|Ascc3   |ENSMUSG00000038774 |   1.2473903| 4.150904| 0.7640501|  0.1216579| 0.9998133|10  |   50592669|  50851202|Females |
|Ascc3   |ENSMUSG00000038774 | -10.0214742| 4.441492| 0.0250024|  0.0748261| 0.0660853|10  |   50592669|  50851202|HFD     |
|Ascc3   |ENSMUSG00000038774 |   1.2473903| 4.150904| 0.7640501|  0.1216579| 0.9998133|10  |   50592669|  50851202|Males   |
|Ascc3   |ENSMUSG00000038774 |   1.2473903| 4.150904| 0.7640501|  0.1216579| 0.9998133|10  |   50592669|  50851202|NCD     |
|Bdh1    |ENSMUSG00000046598 |   1.1657501| 3.023389| 0.6999846| -0.0018146| 0.7771737|16  |   31422280|  31458901|All     |
|Bdh1    |ENSMUSG00000046598 |  -0.8035169| 3.785950| 0.8321042|  0.1214909| 0.9998133|16  |   31422280|  31458901|Females |
|Bdh1    |ENSMUSG00000046598 | -11.6191711| 4.390242| 0.0087000|  0.0823910| 0.0457632|16  |   31422280|  31458901|HFD     |
|Bdh1    |ENSMUSG00000046598 |  -0.8035169| 3.785950| 0.8321042|  0.1214909| 0.9998133|16  |   31422280|  31458901|Males   |
|Bdh1    |ENSMUSG00000046598 |  -0.8035169| 3.785950| 0.8321042|  0.1214909| 0.9998133|16  |   31422280|  31458901|NCD     |
|Cdkal1  |ENSMUSG00000006191 | -11.0765308| 3.200005| 0.0005865|  0.0228312| 0.0032304|13  |   29191746|  29855674|All     |
|Cdkal1  |ENSMUSG00000006191 |   2.0826909| 3.207723| 0.5167874|  0.1228782| 0.9998133|13  |   29191746|  29855674|Females |
|Cdkal1  |ENSMUSG00000006191 | -10.6646608| 4.151583| 0.0108455|  0.0807977| 0.0482649|13  |   29191746|  29855674|HFD     |
|Cdkal1  |ENSMUSG00000006191 |   2.0826909| 3.207723| 0.5167874|  0.1228782| 0.9998133|13  |   29191746|  29855674|Males   |
|Cdkal1  |ENSMUSG00000006191 |   2.0826909| 3.207723| 0.5167874|  0.1228782| 0.9998133|13  |   29191746|  29855674|NCD     |
|Chrna2  |ENSMUSG00000022041 |   9.4887506| 1.033720| 0.0000000|  0.1504871| 0.0000000|14  |   66140960|  66152948|All     |
|Chrna2  |ENSMUSG00000022041 |   2.6658706| 1.426189| 0.0628194|  0.1340375| 0.9998133|14  |   66140960|  66152948|Females |
|Chrna2  |ENSMUSG00000022041 |   0.9268317| 1.808669| 0.6088423|  0.0551697| 0.6712439|14  |   66140960|  66152948|HFD     |
|Chrna2  |ENSMUSG00000022041 |   2.6658706| 1.426189| 0.0628194|  0.1340375| 0.9998133|14  |   66140960|  66152948|Males   |
|Chrna2  |ENSMUSG00000022041 |   2.6658706| 1.426189| 0.0628194|  0.1340375| 0.9998133|14  |   66140960|  66152948|NCD     |
|Clstn3  |ENSMUSG00000008153 |   8.4417805| 1.086837| 0.0000000|  0.1120864| 0.0000000|6   |  124430759| 124464794|All     |
|Clstn3  |ENSMUSG00000008153 |   2.0198754| 1.443326| 0.1629775|  0.1284961| 0.9998133|6   |  124430759| 124464794|Females |
|Clstn3  |ENSMUSG00000008153 |   1.8039296| 1.597981| 0.2601399|  0.0593575| 0.3295498|6   |  124430759| 124464794|HFD     |
|Clstn3  |ENSMUSG00000008153 |   2.0198754| 1.443326| 0.1629775|  0.1284961| 0.9998133|6   |  124430759| 124464794|Males   |
|Clstn3  |ENSMUSG00000008153 |   2.0198754| 1.443326| 0.1629775|  0.1284961| 0.9998133|6   |  124430759| 124464794|NCD     |
|Hmgcl   |ENSMUSG00000028672 |   6.6256980| 3.197128| 0.0387748|  0.0069614| 0.0849394|4   |  135946448| 135962617|All     |
|Hmgcl   |ENSMUSG00000028672 |   4.4436692| 3.463927| 0.2007955|  0.1273586| 0.9998133|4   |  135946448| 135962617|Females |
|Hmgcl   |ENSMUSG00000028672 |  -7.3761984| 4.414412| 0.0961124|  0.0655700| 0.1503337|4   |  135946448| 135962617|HFD     |
|Hmgcl   |ENSMUSG00000028672 |   4.4436692| 3.463927| 0.2007955|  0.1273586| 0.9998133|4   |  135946448| 135962617|Males   |
|Hmgcl   |ENSMUSG00000028672 |   4.4436692| 3.463927| 0.2007955|  0.1273586| 0.9998133|4   |  135946448| 135962617|NCD     |
|Hmgcr   |ENSMUSG00000021670 |  -9.3235352| 1.478871| 0.0000000|  0.0761610| 0.0000001|13  |   96648967|  96670936|All     |
|Hmgcr   |ENSMUSG00000021670 |  -1.9593126| 1.715971| 0.2546806|  0.1261116| 0.9998133|13  |   96648967|  96670936|Females |
|Hmgcr   |ENSMUSG00000021670 |  -1.6616427| 2.234514| 0.4578719|  0.0563755| 0.5284518|13  |   96648967|  96670936|HFD     |
|Hmgcr   |ENSMUSG00000021670 |  -1.9593126| 1.715971| 0.2546806|  0.1261116| 0.9998133|13  |   96648967|  96670936|Males   |
|Hmgcr   |ENSMUSG00000021670 |  -1.9593126| 1.715971| 0.2546806|  0.1261116| 0.9998133|13  |   96648967|  96670936|NCD     |
|Hmgcs2  |ENSMUSG00000027875 |  -3.0947073| 3.298506| 0.3486169| -0.0002549| 0.4623869|3   |   98280435|  98310738|All     |
|Hmgcs2  |ENSMUSG00000027875 |  -1.6768248| 3.617447| 0.6434023|  0.1221172| 0.9998133|3   |   98280435|  98310738|Females |
|Hmgcs2  |ENSMUSG00000027875 |  -9.9989115| 4.343364| 0.0222350|  0.0756572| 0.0626043|3   |   98280435|  98310738|HFD     |
|Hmgcs2  |ENSMUSG00000027875 |  -1.6768248| 3.617447| 0.6434023|  0.1221172| 0.9998133|3   |   98280435|  98310738|Males   |
|Hmgcs2  |ENSMUSG00000027875 |  -1.6768248| 3.617447| 0.6434023|  0.1221172| 0.9998133|3   |   98280435|  98310738|NCD     |
|Hnf1a   |ENSMUSG00000029556 | -16.3148558| 3.917068| 0.0000371|  0.0336134| 0.0003847|5   |  114948980| 114971094|All     |
|Hnf1a   |ENSMUSG00000029556 |   4.5251813| 3.936726| 0.2515125|  0.1261758| 0.9998133|5   |  114948980| 114971094|Females |
|Hnf1a   |ENSMUSG00000029556 |  -9.0876722| 5.508187| 0.1003570|  0.0652851| 0.1551037|5   |  114948980| 114971094|HFD     |
|Hnf1a   |ENSMUSG00000029556 |   4.5251813| 3.936726| 0.2515125|  0.1261758| 0.9998133|5   |  114948980| 114971094|Males   |
|Hnf1a   |ENSMUSG00000029556 |   4.5251813| 3.936726| 0.2515125|  0.1261758| 0.9998133|5   |  114948980| 114971094|NCD     |
|Kif13b  |ENSMUSG00000060012 | -15.0990849| 2.805736| 0.0000001|  0.0561502| 0.0000043|14  |   64652531|  64806296|All     |
|Kif13b  |ENSMUSG00000060012 |  -5.2070985| 3.479772| 0.1358767|  0.1295144| 0.9998133|14  |   64652531|  64806296|Females |
|Kif13b  |ENSMUSG00000060012 |  -7.4887100| 3.647610| 0.0412147|  0.0713207| 0.0856758|14  |   64652531|  64806296|HFD     |
|Kif13b  |ENSMUSG00000060012 |  -5.2070985| 3.479772| 0.1358767|  0.1295144| 0.9998133|14  |   64652531|  64806296|Males   |
|Kif13b  |ENSMUSG00000060012 |  -5.2070985| 3.479772| 0.1358767|  0.1295144| 0.9998133|14  |   64652531|  64806296|NCD     |
|Ldlr    |ENSMUSG00000032193 | -13.3166539| 2.803864| 0.0000027|  0.0438540| 0.0000510|9   |   21723576|  21749916|All     |
|Ldlr    |ENSMUSG00000032193 |  -3.5246181| 3.454263| 0.3085901|  0.1251517| 0.9998133|9   |   21723576|  21749916|Females |
|Ldlr    |ENSMUSG00000032193 |  -3.1375886| 3.696645| 0.3969054|  0.0570692| 0.4683115|9   |   21723576|  21749916|HFD     |
|Ldlr    |ENSMUSG00000032193 |  -3.5246181| 3.454263| 0.3085901|  0.1251517| 0.9998133|9   |   21723576|  21749916|Males   |
|Ldlr    |ENSMUSG00000032193 |  -3.5246181| 3.454263| 0.3085901|  0.1251517| 0.9998133|9   |   21723576|  21749916|NCD     |
|Lpin1   |ENSMUSG00000020593 |  -1.1472997| 1.534029| 0.4548955| -0.0009384| 0.5656342|12  |   16535669|  16589770|All     |
|Lpin1   |ENSMUSG00000020593 |  -2.4975215| 1.467126| 0.0900005|  0.1318947| 0.9998133|12  |   16535669|  16589770|Females |
|Lpin1   |ENSMUSG00000020593 |   0.2900484| 2.044354| 0.8873024|  0.0541606| 0.9099546|12  |   16535669|  16589770|HFD     |
|Lpin1   |ENSMUSG00000020593 |  -2.4975215| 1.467126| 0.0900005|  0.1318947| 0.9998133|12  |   16535669|  16589770|Males   |
|Lpin1   |ENSMUSG00000020593 |  -2.4975215| 1.467126| 0.0900005|  0.1318947| 0.9998133|12  |   16535669|  16589770|NCD     |
|Mmab    |ENSMUSG00000029575 | -18.5419930| 2.312221| 0.0000000|  0.1187055| 0.0000000|5   |  114431034| 114444060|All     |
|Mmab    |ENSMUSG00000029575 |   0.8202457| 2.765015| 0.7669913|  0.1216494| 0.9998133|5   |  114431034| 114444060|Females |
|Mmab    |ENSMUSG00000029575 |  -8.6809209| 3.629687| 0.0175885|  0.0773264| 0.0567946|5   |  114431034| 114444060|HFD     |
|Mmab    |ENSMUSG00000029575 |   0.8202457| 2.765015| 0.7669913|  0.1216494| 0.9998133|5   |  114431034| 114444060|Males   |
|Mmab    |ENSMUSG00000029575 |   0.8202457| 2.765015| 0.7669913|  0.1216494| 0.9998133|5   |  114431034| 114444060|NCD     |
|Mvk     |ENSMUSG00000041939 | -16.1407894| 1.847984| 0.0000000|  0.1380694| 0.0000000|5   |  114444269| 114460591|All     |
|Mvk     |ENSMUSG00000041939 |   2.8851504| 2.371010| 0.2248687|  0.1267575| 0.9998133|5   |  114444269| 114460591|Females |
|Mvk     |ENSMUSG00000041939 |  -5.3849461| 3.270013| 0.1009899|  0.0652438| 0.1558069|5   |  114444269| 114460591|HFD     |
|Mvk     |ENSMUSG00000041939 |   2.8851504| 2.371010| 0.2248687|  0.1267575| 0.9998133|5   |  114444269| 114460591|Males   |
|Mvk     |ENSMUSG00000041939 |   2.8851504| 2.371010| 0.2248687|  0.1267575| 0.9998133|5   |  114444269| 114460591|NCD     |
|Slc16a1 |ENSMUSG00000032902 |  -7.6480159| 2.500941| 0.0023552|  0.0174593| 0.0095752|3   |  104638668| 104658462|All     |
|Slc16a1 |ENSMUSG00000032902 |   3.8790093| 2.424642| 0.1109634|  0.1306733| 0.9998133|3   |  104638668| 104658462|Females |
|Slc16a1 |ENSMUSG00000032902 |  -9.0781974| 3.304002| 0.0064848|  0.0845235| 0.0439051|3   |  104638668| 104658462|HFD     |
|Slc16a1 |ENSMUSG00000032902 |   3.8790093| 2.424642| 0.1109634|  0.1306733| 0.9998133|3   |  104638668| 104658462|Males   |
|Slc16a1 |ENSMUSG00000032902 |   3.8790093| 2.424642| 0.1109634|  0.1306733| 0.9998133|3   |  104638668| 104658462|NCD     |
|Slc16a6 |ENSMUSG00000041920 |  -4.8413216| 2.262396| 0.0328775|  0.0075578| 0.0745890|11  |  109450855| 109473598|All     |
|Slc16a6 |ENSMUSG00000041920 |   1.5117865| 2.059688| 0.4636784|  0.1233091| 0.9998133|11  |  109450855| 109473598|Females |
|Slc16a6 |ENSMUSG00000041920 |  -5.0582510| 3.260100| 0.1221598|  0.0640031| 0.1805235|11  |  109450855| 109473598|HFD     |
|Slc16a6 |ENSMUSG00000041920 |   1.5117865| 2.059688| 0.4636784|  0.1233091| 0.9998133|11  |  109450855| 109473598|Males   |
|Slc16a6 |ENSMUSG00000041920 |   1.5117865| 2.059688| 0.4636784|  0.1233091| 0.9998133|11  |  109450855| 109473598|NCD     |
|Slc16a7 |ENSMUSG00000020102 |  -3.9297566| 1.804177| 0.0298928|  0.0079036| 0.0693780|10  |  125226464| 125328963|All     |
|Slc16a7 |ENSMUSG00000020102 |  -3.9519130| 2.548201| 0.1222631|  0.1301155| 0.9998133|10  |  125226464| 125328963|Females |
|Slc16a7 |ENSMUSG00000020102 |  -1.3556743| 3.283676| 0.6801043|  0.0547865| 0.7343205|10  |  125226464| 125328963|HFD     |
|Slc16a7 |ENSMUSG00000020102 |  -3.9519130| 2.548201| 0.1222631|  0.1301155| 0.9998133|10  |  125226464| 125328963|Males   |
|Slc16a7 |ENSMUSG00000020102 |  -3.9519130| 2.548201| 0.1222631|  0.1301155| 0.9998133|10  |  125226464| 125328963|NCD     |
|Srebf1  |ENSMUSG00000020538 |   4.5375786| 2.110031| 0.0320278|  0.0076528| 0.0731406|11  |   60199089|  60222581|All     |
|Srebf1  |ENSMUSG00000020538 |   1.0428632| 2.269138| 0.6462331|  0.1221037| 0.9998133|11  |   60199089|  60222581|Females |
|Srebf1  |ENSMUSG00000020538 |  -0.4879524| 3.077482| 0.8741598|  0.0541815| 0.8998479|11  |   60199089|  60222581|HFD     |
|Srebf1  |ENSMUSG00000020538 |   1.0428632| 2.269138| 0.6462331|  0.1221037| 0.9998133|11  |   60199089|  60222581|Males   |
|Srebf1  |ENSMUSG00000020538 |   1.0428632| 2.269138| 0.6462331|  0.1221037| 0.9998133|11  |   60199089|  60222581|NCD     |
:::
:::


### Cholesterol Data Heatmap


::: {.cell}

```{.r .cell-code}
library(pheatmap)
cholesterol.genes <- c('Hmgcl','Hmgcs2','Hmgcr','Bdh1','Slc16a1','Slc16a3','Slc16a6','Slc16a7','Slc16a8','Ldlr','Mvk')

heatmap.data <- filter(summary.data, symbol %in% cholesterol.genes) %>%
  select(beta,model,symbol) %>% 
  pivot_wider(names_from=model,values_from=beta) %>%
  column_to_rownames('symbol') %>%
  select(order(column.order)) %>%
  as.matrix

#to make coloring symmetrical
max_value <- max(abs(heatmap.data))
breaks <- seq(-max_value, max_value, length.out = 101)

pheatmap(heatmap.data,
         main="Effect Estimates",
         fontsize=16,
         breaks = breaks,
         color = colorRampPalette(c("#00274C", "white", "#FFCB05"))(100),
         #cellwidth=20,
         #cellheight=20,
         cluster_rows=T,
         cluster_cols=F,
         show_rownames=T,
         show_colnames=T)
```

::: {.cell-output-display}
![](liver-mRNA-cholesterol-analyses_files/figure-html/cholestrol-heatmap-1.png){width=672}
:::
:::


# Session Information


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}
```
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

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
 [1] pheatmap_1.0.12 knitr_1.49      janitor_2.2.1   lubridate_1.9.4
 [5] forcats_1.0.0   stringr_1.5.1   purrr_1.0.2     tibble_3.2.1   
 [9] ggplot2_3.5.1   tidyverse_2.0.0 broom_1.0.7     tidyr_1.3.1    
[13] dplyr_1.1.4     readr_2.1.5    

loaded via a namespace (and not attached):
 [1] generics_0.1.3     stringi_1.8.4      lattice_0.22-6     hms_1.1.3         
 [5] digest_0.6.37      magrittr_2.0.3     RColorBrewer_1.1-3 evaluate_1.0.1    
 [9] grid_4.4.0         timechange_0.3.0   fastmap_1.2.0      jsonlite_1.8.9    
[13] Matrix_1.7-0       backports_1.5.0    mgcv_1.9-1         scales_1.3.0      
[17] textshaping_0.4.1  cli_3.6.3          crayon_1.5.3       rlang_1.1.4       
[21] bit64_4.5.2        munsell_0.5.1      splines_4.4.0      withr_3.0.2       
[25] yaml_2.3.9         parallel_4.4.0     tools_4.4.0        tzdb_0.4.0        
[29] colorspace_2.1-1   vctrs_0.6.5        R6_2.5.1           lifecycle_1.0.4   
[33] snakecase_0.11.1   bit_4.5.0.1        htmlwidgets_1.6.4  vroom_1.6.5       
[37] ragg_1.3.3         pkgconfig_2.0.3    pillar_1.10.1      gtable_0.3.6      
[41] glue_1.8.0         systemfonts_1.1.0  xfun_0.50          tidyselect_1.2.1  
[45] rstudioapi_0.17.1  farver_2.1.2       htmltools_0.5.8.1  nlme_3.1-164      
[49] rmarkdown_2.29     labeling_0.4.3     compiler_4.4.0    
```
:::
:::
