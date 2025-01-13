---
title: "Analysis of Diversity Outbred Strain Cholesterol Levels Separating by Diet"
author: "Dave Bridges"
date: "July 7, 2020"
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

# Experimental Details

This analyses the data analysed via GEMMA and provided in the various output folders.

# Data Entry


``` r
# additive lmm data
lmm.additive.filename <- 'output/cholesterol_all.assoc.txt'
lmm.additive.data <- read_tsv(lmm.additive.filename) %>%
  separate(rs, sep="_", into=c("chromosome","position","allele","alt"),remove=FALSE) %>%
  mutate(chromosome=factor(chromosome, c(1:19,"X"))) %>%
  mutate(position=as.integer(position))

lmm.ncd.filename <- 'output/cholesterol_ncd.assoc.txt'
lmm.ncd.results <- read_tsv(lmm.ncd.filename) %>%
  separate(rs, sep="_", into=c("chromosome","position","allele","alt"),remove=FALSE) %>%
  mutate(chromosome=factor(chromosome, c(1:19,"X"))) %>%
  mutate(position=as.integer(position))

lmm.hfd.filename <- 'output/cholesterol_hfd.assoc.txt'
lmm.hfd.results <- read_tsv(lmm.hfd.filename) %>%
  separate(rs, sep="_", into=c("chromosome","position","allele","alt"),remove=FALSE) %>%
  mutate(chromosome=factor(chromosome, c(1:19,"X"))) %>%
  mutate(position=as.integer(position))
```


# Additive Models

Ran the additive models using GEMMA, first using intercepts and additiv ecovariates for diet and sex

$SNP = beta_1 SNP + \beta_2 Diet + \beta_3 Sex + \mu +\epsilon$

Where $$\epsilon$$ are the residuals and $$\mu$$ is the relationship matrix of the strains as defined by

## LMM Analysis for Additive Models


``` r
library(qqman)
qq(lmm.additive.data$p_wald)
suggestive.pval <- 1E-5
genome.pval <- 5E-8

lmm.additive.data %>%
  arrange(p_wald) %>% 
  filter(p_wald<genome.pval) %>%
  kable(caption="Genome-wide significant associations from mixed linear models for cholesterol in additive model") 
```



Table: Genome-wide significant associations from mixed linear models for cholesterol in additive model

|chr |rs              |chromosome |  position|allele |alt |        ps| n_miss|allele1 |allele0 |    af| beta|   se| logl_H1| l_remle| p_wald|
|:---|:---------------|:----------|---------:|:------|:---|---------:|------:|:-------|:-------|-----:|----:|----:|-------:|-------:|------:|
|1   |1_171425406_H_C |1          | 171425406|H      |C   | 171425406|      0|C       |H       | 0.094| 14.6| 2.62|   -3815|    2.22|      0|
|1   |1_171431917_H_C |1          | 171431917|H      |C   | 171431917|      0|C       |H       | 0.094| 14.6| 2.62|   -3815|    2.22|      0|
|1   |1_171295194_H_C |1          | 171295194|H      |C   | 171295194|      0|C       |H       | 0.095| 14.5| 2.61|   -3815|    2.21|      0|
|1   |1_171418895_H_C |1          | 171418895|H      |C   | 171418895|      0|C       |H       | 0.095| 14.5| 2.62|   -3815|    2.22|      0|

``` r
lmm.additive.data %>%
  arrange(p_wald) %>% 
  filter(p_wald<suggestive.pval) %>%
  mutate(position.start = substr(as.character(position), 1,2)) %>%
  group_by(chromosome,position.start) %>%
  summarize_all(.funs=first) %>%
  select(-position.start,-chr,-ps) %>%
  kable(caption="Suggestive genome-wide significant associations from mixed linear models for cholesterol in additive model, clumped by first two digits of the position") 
```



Table: Suggestive genome-wide significant associations from mixed linear models for cholesterol in additive model, clumped by first two digits of the position

|chromosome |rs              |  position|allele |alt | n_miss|allele1 |allele0 |    af|  beta|   se| logl_H1| l_remle| p_wald|
|:----------|:---------------|---------:|:------|:---|------:|:-------|:-------|-----:|-----:|----:|-------:|-------:|------:|
|1          |1_138064771_D_H | 138064771|D      |H   |      0|H       |D       | 0.160| 10.41| 2.14|   -3818|    1.75|      0|
|1          |1_143704476_D_H | 143704476|D      |H   |      0|H       |D       | 0.161|  9.94| 2.15|   -3820|    1.79|      0|
|1          |1_155276305_H_H | 155276305|H      |H   |      0|H       |H       | 0.199|  9.83| 1.98|   -3818|    1.64|      0|
|1          |1_169972414_F_C | 169972414|F      |C   |      0|C       |F       | 0.084| 14.85| 2.72|   -3815|    2.23|      0|
|1          |1_171425406_H_C | 171425406|H      |C   |      0|C       |H       | 0.094| 14.59| 2.62|   -3815|    2.22|      0|
|5          |5_123629774_B_E | 123629774|B      |E   |      0|E       |B       | 0.110| 11.66| 2.45|   -3819|    2.17|      0|

``` r
lmm.additive.data %>%
  arrange(p_wald) %>% 
  filter(p_wald<1E-4) %>%
  mutate(position.start = substr(as.character(position), 1,2)) %>%
  group_by(chromosome,position.start) %>%
  summarize_all(.funs=first) %>%
  select(-position.start,-chr,-ps) -> additive.snp.summary

additive.snp.summary %>%
  kable(caption="Relaxed suggestive genome-wide significant associations from mixed linear models for cholesterol in additive model, clumped by first two digits of the position") 
```



Table: Relaxed suggestive genome-wide significant associations from mixed linear models for cholesterol in additive model, clumped by first two digits of the position

|chromosome |rs              |  position|allele |alt | n_miss|allele1 |allele0 |    af|  beta|   se| logl_H1| l_remle| p_wald|
|:----------|:---------------|---------:|:------|:---|------:|:-------|:-------|-----:|-----:|----:|-------:|-------:|------:|
|1          |1_138064771_D_H | 138064771|D      |H   |      0|H       |D       | 0.160| 10.41| 2.14|   -3818|    1.75|      0|
|1          |1_143704476_D_H | 143704476|D      |H   |      0|H       |D       | 0.161|  9.94| 2.15|   -3820|    1.79|      0|
|1          |1_155276305_H_H | 155276305|H      |H   |      0|H       |H       | 0.199|  9.83| 1.98|   -3818|    1.64|      0|
|1          |1_169972414_F_C | 169972414|F      |C   |      0|C       |F       | 0.084| 14.85| 2.72|   -3815|    2.23|      0|
|1          |1_171425406_H_C | 171425406|H      |C   |      0|C       |H       | 0.094| 14.59| 2.62|   -3815|    2.22|      0|
|3          |3_148547289_A_H | 148547289|A      |H   |      0|H       |A       | 0.109| 10.00| 2.48|   -3822|    1.99|      0|
|5          |5_123629774_B_E | 123629774|B      |E   |      0|E       |B       | 0.110| 11.66| 2.45|   -3819|    2.17|      0|
|13         |13_29976363_E_C |  29976363|E      |C   |      0|C       |E       | 0.156|  9.25| 2.19|   -3821|    2.10|      0|
|13         |13_30180778_E_C |  30180778|E      |C   |      0|C       |E       | 0.156|  9.27| 2.19|   -3821|    2.10|      0|

``` r
library(ggmanh)
```

![](figures-diet/lmm-additive-1.png)<!-- -->

``` r
manhattan_plot(x = lmm.additive.data, pval.colname = "p_wald", chr.colname = "chromosome", pos.colname = "position", plot.title = "DO Mice on HFHS Diet (Additive Model)", y.label = "LOD Score")
```

![](figures-diet/lmm-additive-2.png)<!-- -->

``` r
#library(GWASTools)
#with(lmm.additive.data,manhattanPlot(p=p_wald,
                               #chromosome=chromosome,
                               #signif=suggestive.pval))
```

### Chromosome 1 Peak


``` r
snp.pos <- 171425406
library(forcats)
peak1.b <- filter(lmm.additive.data,
                  chromosome==1,
                  position>snp.pos-20000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak1.b,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 1") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=8))
```

![](figures-diet/cholesterol-chr1-160000000-185000000-all-1.png)<!-- -->

``` r
snp.pos <- 80237174
library(forcats)
peak1 <- filter(lmm.additive.data,
                  chromosome==1,
                  position>snp.pos-35000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak1,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 1") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=4))
```

![](figures-diet/cholesterol-chr1-60000000-100000000-all-1.png)<!-- -->

``` r
snp.pos <- 143704476
library(forcats)
peak1 <- filter(lmm.additive.data,
                  chromosome==1,
                  position>snp.pos-20000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak1,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 1") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=4))
```

![](figures-diet/cholesterol-chr1-130000000-160000000-all-1.png)<!-- -->



``` r
snp.pos <- 148547289 
library(forcats)
peak3 <- filter(lmm.additive.data,
                  chromosome==3,
                  position>snp.pos-20000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak3,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 3") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=4))
```

![](figures-diet/cholesterol-chr3-130000000-160000000-all-1.png)<!-- -->




``` r
snp.pos <- 16902424
library(forcats)
peak2 <- filter(lmm.additive.data,
                  chromosome==2,
                  position>snp.pos-20000000,
                  position<snp.pos+10000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak2,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 2") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=5))
```

![](figures-diet/cholesterol-chr2--all-1.png)<!-- -->


``` r
snp.pos <-55298460

library(forcats)
peak5 <- filter(lmm.additive.data,
                  chromosome==5,
                  position>snp.pos-20000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak5,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 5") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=4))
```

![](figures-diet/cholesterol-chr5-60000000-100000000-all-1.png)<!-- -->



### Chromosome 5 Peak


``` r
snp.pos <- 123629774
peak5 <- filter(lmm.additive.data,
                  chromosome==5,
                  position>snp.pos-10000000,
                  position<snp.pos+10000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak5,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 5") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=8))
```

![](figures-diet/chr5-120000000-130000000-1.png)<!-- -->
##Chromosome 10

``` r
snp.pos <- 54344829
peak10 <- filter(lmm.additive.data,
                  chromosome==10,
                  position>snp.pos-20000000,
                  position<snp.pos+30000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak10,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 10") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=8))
```

![](figures-diet/chr10-50000000-80000000-1.png)<!-- -->
##Chromosome 10

``` r
snp.pos <- 95382136
library(forcats)
peak10.b <- filter(lmm.additive.data,
                  chromosome==10,
                  position>snp.pos-20000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak10.b,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 10") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=8))
```

![](figures-diet/cholesterol-chr10-85000000-110000000-all-1.png)<!-- -->























``` r
snp.pos <- 91658214
peak11 <- filter(lmm.additive.data,
                  chromosome==11,
                  position>snp.pos-20000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak11,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 11") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=8))
```

![](figures-diet/chr11-80000000-110000000-1.png)<!-- -->




``` r
snp.pos <-19934709
peak14 <- filter(lmm.additive.data,
                  chromosome==14,
                  position>snp.pos-40000000,
                  position<snp.pos+10000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                       "WSB"="H"))
# 
# ggplot(data=peak14,
#        aes(x=position,
#        y=beta,
#        col=alt,
#        group=alt)) +
#   geom_line() +
#   geom_hline(yintercept=0,lty=2) +
#   labs(title="Strain Specific Effects",
#        y="Cholesterol Effect Size (mg/dL)",
#        x="Position on Chromosome 14") +
#   scale_color_discrete(name="") +
#   guides(col=guide_legend(ncol=2)) +
#   theme_classic(base_size=12) +
#   theme(legend.position=c(0.18,0.9),
#         legend.text=element_text(size=5))
```

``` r
library(readr)
#figures makde will go to directory called figures, will make them as both png and pdf files 
opts_chunk$set(fig.path='figures/',
               echo=TRUE, warning=FALSE, message=FALSE,dev=c('png','pdf'))
options(scipen = 2, digits = 3)
snp.pos <- 4.5e07
library(forcats)
peak9 <- filter(lmm.additive.data,
                  chromosome==9,
                  position>snp.pos-20000000,
                  position<snp.pos+20000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak9,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 9") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=5))
```

![](figures-diet/cholesterol-chr9--all-1.png)<!-- -->


##Chromosome 17

``` r
snp.pos <- 80080939
peak17 <- filter(lmm.additive.data,
                  chromosome==17,
                  position>snp.pos-10000000,
                  position<snp.pos+10000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

 ggplot(data=peak17,
        aes(x=position,        
        y=beta,
        col=alt,
        group=alt)) +
   geom_line() +
   geom_hline(yintercept=0,lty=2) +
   labs(title="Strain Specific Effects",
        y="Cholesterol Effect Size (mg/dL)",
        x="Position on Chromosome 17") +
   scale_color_discrete(name="") +
   guides(col=guide_legend(ncol=2)) +
   theme_classic(base_size=12) +
   theme(legend.position=c(0.18,0.9),
         legend.text=element_text(size=5))
```

![](figures/chr17-75000000-90000000-1.png)<!-- -->





``` r
library(readr)
#figures makde will go to directory called figures, will make them as both png and pdf files 
opts_chunk$set(fig.path='figures/',
               echo=TRUE, warning=FALSE, message=FALSE,dev=c('png','pdf'))
options(scipen = 2, digits = 3)
snp.pos <- 2e07
library(forcats)
peak12 <- filter(lmm.additive.data,
                  chromosome==12,
                  position>snp.pos-1e07,
                  position<snp.pos+1e07) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H"))

ggplot(data=peak12,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size (mg/dL)",
       x="Position on Chromosome 12") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.18,0.9),
        legend.text=element_text(size=5))
```

![](figures/cholesterol-chr12--all-1.png)<!-- -->

## SNP Analysis

### Linear Mixed Model SNP Analysis for Chow


``` r
qq(lmm.ncd.results$p_wald)
```

![](figures/chol-snp-analysis-chow-1.png)<!-- -->

``` r
lmm.ncd.results %>%
  arrange(p_wald) %>% 
  filter(p_wald<genome.pval) %>%
  kable(caption="Genome-wide significant associations from mixed linear models for cholesterol on NCD") 
```



Table: Genome-wide significant associations from mixed linear models for cholesterol on NCD

|chr |rs |chromosome | position|allele |alt | ps| n_miss|allele1 |allele0 | af| beta| se| logl_H1| l_remle| p_wald|
|:---|:--|:----------|--------:|:------|:---|--:|------:|:-------|:-------|--:|----:|--:|-------:|-------:|------:|

``` r
lmm.ncd.results %>%
   arrange(p_wald) %>% 
   filter(p_wald<suggestive.pval) %>%
  mutate(position.start = substr(as.character(position), 1,2)) %>%
  group_by(chromosome,position.start) %>%
  summarize_all(.funs=first) %>%
  select(-position.start,-chr,-ps) -> ncd.snp.summary
  
lmm.ncd.results %>%
   arrange(p_wald) %>% 
   filter(p_wald<1E-3) %>%
   mutate(position.start = substr(as.character(position), 1,2)) %>%
  group_by(chromosome,position.start) %>%
  summarize_all(.funs=first) %>%
  select(-position.start,-chr,-ps)  %>% 
  kable(caption="Relaxed suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position") 
```



Table: Relaxed suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position

|chromosome |rs               |  position|allele |alt | n_miss|allele1 |allele0 |    af|   beta|   se| logl_H1| l_remle| p_wald|
|:----------|:----------------|---------:|:------|:---|------:|:-------|:-------|-----:|------:|----:|-------:|-------:|------:|
|1          |1_117784609_B_C  | 117784609|B      |C   |      0|C       |B       | 0.119|   9.25| 2.62|   -1966|    3.30|  0.000|
|1          |1_136644599_C_A  | 136644599|C      |A   |      0|A       |C       | 0.151|  -9.19| 2.42|   -1966|    2.51|  0.000|
|1          |1_152428168_H_A  | 152428168|H      |A   |      0|A       |H       | 0.139|  -8.99| 2.45|   -1966|    3.07|  0.000|
|1          |1_169972414_F_C  | 169972414|F      |C   |      0|C       |F       | 0.087|  14.49| 2.93|   -1960|    3.55|  0.000|
|1          |1_170324641_H_C  | 170324641|H      |C   |      0|C       |H       | 0.087|  14.52| 2.93|   -1960|    3.55|  0.000|
|1          |1_183576795_H_C  | 183576795|H      |C   |      0|C       |H       | 0.099|   9.30| 2.74|   -1966|    3.68|  0.001|
|2          |2_14535771_H_G   |  14535771|H      |G   |      0|G       |H       | 0.132|   8.58| 2.51|   -1967|    3.06|  0.001|
|2          |2_15491718_H_G   |  15491718|H      |G   |      0|G       |H       | 0.111|  10.14| 2.64|   -1965|    2.80|  0.000|
|2          |2_16663582_H_G   |  16663582|H      |G   |      0|G       |H       | 0.110|  10.47| 2.65|   -1965|    2.83|  0.000|
|2          |2_17009990_H_G   |  17009990|H      |G   |      0|G       |H       | 0.115|   9.50| 2.60|   -1966|    2.80|  0.000|
|2          |2_18921607_H_G   |  18921607|H      |G   |      0|G       |H       | 0.115|   9.31| 2.56|   -1966|    2.92|  0.000|
|2          |2_19865037_H_G   |  19865037|H      |G   |      0|G       |H       | 0.114|   9.19| 2.62|   -1966|    3.12|  0.000|
|2          |2_20154521_H_G   |  20154521|H      |G   |      0|G       |H       | 0.109|   8.70| 2.58|   -1967|    3.20|  0.001|
|2          |2_43690499_H_G   |  43690499|H      |G   |      0|G       |H       | 0.113|   9.66| 2.60|   -1965|    3.70|  0.000|
|2          |2_44240389_H_G   |  44240389|H      |G   |      0|G       |H       | 0.119|   8.89| 2.49|   -1966|    3.73|  0.000|
|2          |2_45242400_H_G   |  45242400|H      |G   |      0|G       |H       | 0.118|   9.18| 2.58|   -1966|    3.70|  0.000|
|2          |2_47871335_G_G   |  47871335|G      |G   |      0|G       |G       | 0.136|   8.04| 2.39|   -1967|    3.56|  0.001|
|2          |2_48818189_G_G   |  48818189|G      |G   |      0|G       |G       | 0.141|   8.53| 2.35|   -1966|    3.54|  0.000|
|2          |2_49650224_G_G   |  49650224|G      |G   |      0|G       |G       | 0.140|   9.14| 2.36|   -1965|    3.61|  0.000|
|3          |3_156886929_F_B  | 156886929|F      |B   |      0|B       |F       | 0.105|   9.06| 2.71|   -1967|    4.00|  0.001|
|4          |4_54753377_C_G   |  54753377|C      |G   |      0|G       |C       | 0.122|   8.77| 2.56|   -1966|    3.30|  0.001|
|5          |5_128399127_D_E  | 128399127|D      |E   |      0|E       |D       | 0.118|   9.15| 2.53|   -1966|    3.74|  0.000|
|5          |5_133357791_C_F  | 133357791|C      |F   |      0|F       |C       | 0.096|   9.99| 2.93|   -1966|    3.48|  0.001|
|6          |6_126771534_D_E  | 126771534|D      |E   |      0|E       |D       | 0.065|  12.67| 3.27|   -1965|    3.67|  0.000|
|6          |6_41772867_G_G   |  41772867|G      |G   |      0|G       |G       | 0.119|  -9.01| 2.57|   -1966|    3.60|  0.001|
|6          |6_42232037_G_G   |  42232037|G      |G   |      0|G       |G       | 0.118|  -9.02| 2.56|   -1966|    3.62|  0.000|
|10         |10_54747088_H_B  |  54747088|H      |B   |      0|B       |H       | 0.065| -10.96| 3.05|   -1966|    3.25|  0.000|
|10         |10_55538600_H_B  |  55538600|H      |B   |      0|B       |H       | 0.063| -11.75| 3.13|   -1965|    3.23|  0.000|
|11         |11_121921138_G_G | 121921138|G      |G   |      0|G       |G       | 0.076|  -9.88| 2.97|   -1967|    3.65|  0.001|
|12         |12_84203316_H_C  |  84203316|H      |C   |      0|C       |H       | 0.206|  -7.26| 2.13|   -1966|    3.75|  0.001|
|14         |14_58887010_F_F  |  58887010|F      |F   |      0|F       |F       | 0.108|   9.01| 2.65|   -1966|    3.52|  0.001|
|14         |14_59011591_F_F  |  59011591|F      |F   |      0|F       |F       | 0.107|   9.06| 2.67|   -1966|    3.48|  0.001|
|14         |14_62912822_H_F  |  62912822|H      |F   |      0|F       |H       | 0.096|  10.09| 2.79|   -1966|    3.50|  0.000|
|14         |14_63527073_H_F  |  63527073|H      |F   |      0|F       |H       | 0.101|  10.91| 2.76|   -1964|    3.53|  0.000|
|14         |14_64610241_H_F  |  64610241|H      |F   |      0|F       |H       | 0.102|  11.04| 2.75|   -1964|    3.52|  0.000|
|14         |14_65047597_H_F  |  65047597|H      |F   |      0|F       |H       | 0.103|  11.02| 2.75|   -1964|    3.52|  0.000|
|14         |14_66022446_H_F  |  66022446|H      |F   |      0|F       |H       | 0.103|   9.80| 2.73|   -1966|    3.42|  0.000|
|14         |14_67143276_H_F  |  67143276|H      |F   |      0|F       |H       | 0.090|  10.53| 2.95|   -1966|    3.42|  0.000|
|15         |15_93728085_E_G  |  93728085|E      |G   |      0|G       |E       | 0.058|  12.15| 3.62|   -1967|    3.37|  0.001|
|15         |15_94566593_F_G  |  94566593|F      |G   |      0|G       |F       | 0.062|  12.26| 3.44|   -1966|    3.49|  0.000|
|18         |18_44959139_G_A  |  44959139|G      |A   |      0|A       |G       | 0.148|   8.01| 2.27|   -1966|    3.70|  0.000|
|18         |18_45596845_G_A  |  45596845|G      |A   |      0|A       |G       | 0.145|   8.95| 2.30|   -1965|    3.58|  0.000|
|18         |18_46410922_G_A  |  46410922|G      |A   |      0|A       |G       | 0.136|  11.59| 2.43|   -1961|    3.51|  0.000|
|18         |18_47422415_G_A  |  47422415|G      |A   |      0|A       |G       | 0.135|  10.85| 2.42|   -1962|    3.53|  0.000|
|18         |18_48037460_G_A  |  48037460|G      |A   |      0|A       |G       | 0.135|  10.41| 2.40|   -1963|    3.51|  0.000|
|18         |18_49809849_C_A  |  49809849|C      |A   |      0|A       |C       | 0.146|   8.62| 2.37|   -1966|    3.32|  0.000|
|18         |18_50045754_C_A  |  50045754|C      |A   |      0|A       |C       | 0.144|   8.41| 2.37|   -1966|    3.34|  0.000|

``` r
ncd.snp.summary %>%
  kable(caption="Suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position") 
```



Table: Suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position

|chromosome |rs              |  position|allele |alt | n_miss|allele1 |allele0 |    af| beta|   se| logl_H1| l_remle| p_wald|
|:----------|:---------------|---------:|:------|:---|------:|:-------|:-------|-----:|----:|----:|-------:|-------:|------:|
|1          |1_169972414_F_C | 169972414|F      |C   |      0|C       |F       | 0.087| 14.5| 2.93|   -1960|    3.55|      0|
|1          |1_170324641_H_C | 170324641|H      |C   |      0|C       |H       | 0.087| 14.5| 2.93|   -1960|    3.55|      0|
|18         |18_46410922_G_A |  46410922|G      |A   |      0|A       |G       | 0.136| 11.6| 2.43|   -1961|    3.51|      0|
|18         |18_47422415_G_A |  47422415|G      |A   |      0|A       |G       | 0.135| 10.9| 2.42|   -1962|    3.53|      0|

``` r
manhattan_plot(x = lmm.ncd.results, pval.colname = "p_wald", chr.colname = "chromosome", pos.colname = "position", plot.title = "DO Mice on a NCD Diet", y.label = "LOD Score") -> ncd.manhattan

ncd.manhattan
```

![](figures/chol-snp-analysis-chow-2.png)<!-- -->

### Linear Mixed Model SNP Analysis for HFD


``` r
qq(lmm.hfd.results$p_wald)
```

![](figures/chol-snp-analysis-hfd-1.png)<!-- -->

``` r
lmm.hfd.results %>%
  arrange(p_wald) %>% 
  filter(p_wald<genome.pval) %>%
  kable(caption="Genome-wide significant associations from mixed linear models for cholesterol on HFD") 
```



Table: Genome-wide significant associations from mixed linear models for cholesterol on HFD

|chr |rs |chromosome | position|allele |alt | ps| n_miss|allele1 |allele0 | af| beta| se| logl_H1| l_remle| p_wald|
|:---|:--|:----------|--------:|:------|:---|--:|------:|:-------|:-------|--:|----:|--:|-------:|-------:|------:|

``` r
lmm.hfd.results %>%
   arrange(p_wald) %>% 
   filter(p_wald<suggestive.pval) %>%
   mutate(position.start = substr(as.character(position), 1,2)) %>%
  group_by(chromosome,position.start) %>%
  summarize_all(.funs=first) %>%
  select(-position.start,-chr,-ps) %>%
  kable(caption="Suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position") 
```



Table: Suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position

|chromosome |rs | position|allele |alt | n_miss|allele1 |allele0 | af| beta| se| logl_H1| l_remle| p_wald|
|:----------|:--|--------:|:------|:---|------:|:-------|:-------|--:|----:|--:|-------:|-------:|------:|

``` r
lmm.hfd.results %>%
   arrange(p_wald) %>% 
   filter(p_wald<1E-3) %>%
   mutate(position.start = substr(as.character(position), 1,2)) %>%
  group_by(chromosome,position.start) %>%
  summarize_all(.funs=first) %>%
  select(-position.start,-chr,-ps) ->
  hfhs.snp.summary

hfhs.snp.summary %>%
  kable(caption="Suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position") 
```



Table: Suggestive genome-wide significant associations from mixed linear models for cholesterol on HFD, clumped by first two digits of the position

|chromosome |rs               |  position|allele |alt | n_miss|allele1 |allele0 |    af|  beta|   se| logl_H1| l_remle| p_wald|
|:----------|:----------------|---------:|:------|:---|------:|:-------|:-------|-----:|-----:|----:|-------:|-------:|------:|
|1          |1_137999848_D_H  | 137999848|D      |H   |      0|H       |D       | 0.183|  13.0| 3.39|   -1813|    1.99|  0.000|
|1          |1_145155626_A_H  | 145155626|A      |H   |      0|H       |A       | 0.169|  14.0| 3.45|   -1812|    2.00|  0.000|
|1          |1_155454177_H_H  | 155454177|H      |H   |      0|H       |H       | 0.207|  14.2| 3.19|   -1810|    1.88|  0.000|
|1          |1_171425406_H_C  | 171425406|H      |C   |      0|C       |H       | 0.097|  14.9| 4.42|   -1814|    3.08|  0.001|
|1          |1_83629456_F_B   |  83629456|F      |B   |      0|B       |F       | 0.144| -12.6| 3.74|   -1814|    3.08|  0.001|
|2          |2_44937781_H_E   |  44937781|H      |E   |      0|E       |H       | 0.076| -17.0| 4.97|   -1814|    2.97|  0.001|
|2          |2_45843192_H_E   |  45843192|H      |E   |      0|E       |H       | 0.076| -17.0| 4.95|   -1814|    2.99|  0.001|
|2          |2_46253198_H_E   |  46253198|H      |E   |      0|E       |H       | 0.076| -17.0| 4.95|   -1814|    3.00|  0.001|
|2          |2_47121856_D_E   |  47121856|D      |E   |      0|E       |D       | 0.074| -16.9| 4.97|   -1814|    3.06|  0.001|
|2          |2_5564930_G_E    |   5564930|G      |E   |      0|E       |G       | 0.133| -12.4| 3.73|   -1814|    3.13|  0.001|
|2          |2_5611635_G_E    |   5611635|G      |E   |      0|E       |G       | 0.133| -12.4| 3.73|   -1814|    3.16|  0.001|
|2          |2_5772554_D_E    |   5772554|D      |E   |      0|E       |D       | 0.133| -12.4| 3.70|   -1814|    3.17|  0.001|
|2          |2_6108249_D_E    |   6108249|D      |E   |      0|E       |D       | 0.132| -12.4| 3.70|   -1814|    3.16|  0.001|
|2          |2_6288486_D_E    |   6288486|D      |E   |      0|E       |D       | 0.131| -12.6| 3.70|   -1814|    3.17|  0.001|
|2          |2_6368382_D_E    |   6368382|D      |E   |      0|E       |D       | 0.129| -14.1| 3.72|   -1812|    3.22|  0.000|
|2          |2_6468723_D_E    |   6468723|D      |E   |      0|E       |D       | 0.128| -14.2| 3.73|   -1812|    3.22|  0.000|
|2          |2_6538764_D_E    |   6538764|D      |E   |      0|E       |D       | 0.127| -14.1| 3.73|   -1812|    3.22|  0.000|
|2          |2_6694186_D_E    |   6694186|D      |E   |      0|E       |D       | 0.127| -14.1| 3.73|   -1812|    3.22|  0.000|
|2          |2_6785451_D_E    |   6785451|D      |E   |      0|E       |D       | 0.127| -14.1| 3.74|   -1812|    3.22|  0.000|
|2          |2_6802081_D_E    |   6802081|D      |E   |      0|E       |D       | 0.127| -14.1| 3.74|   -1812|    3.22|  0.000|
|2          |2_7025078_D_E    |   7025078|D      |E   |      0|E       |D       | 0.127| -14.3| 3.75|   -1812|    3.21|  0.000|
|2          |2_7250227_D_E    |   7250227|D      |E   |      0|E       |D       | 0.127| -14.4| 3.76|   -1812|    3.21|  0.000|
|2          |2_7501966_H_E    |   7501966|H      |E   |      0|E       |H       | 0.122| -14.8| 3.82|   -1812|    3.20|  0.000|
|2          |2_7611414_D_E    |   7611414|D      |E   |      0|E       |D       | 0.125| -14.1| 3.80|   -1813|    3.20|  0.000|
|2          |2_7737471_H_E    |   7737471|H      |E   |      0|E       |H       | 0.124| -14.3| 3.81|   -1813|    3.19|  0.000|
|2          |2_7857737_H_E    |   7857737|H      |E   |      0|E       |H       | 0.124| -14.2| 3.81|   -1813|    3.16|  0.000|
|2          |2_8011673_H_E    |   8011673|H      |E   |      0|E       |H       | 0.123| -14.4| 3.82|   -1812|    3.19|  0.000|
|2          |2_8182575_A_E    |   8182575|A      |E   |      0|E       |A       | 0.125| -14.0| 3.78|   -1813|    3.15|  0.000|
|2          |2_8353477_A_E    |   8353477|A      |E   |      0|E       |A       | 0.127| -14.2| 3.77|   -1812|    3.15|  0.000|
|2          |2_8524379_A_E    |   8524379|A      |E   |      0|E       |A       | 0.127| -14.2| 3.76|   -1812|    3.15|  0.000|
|2          |2_8690542_A_E    |   8690542|A      |E   |      0|E       |A       | 0.128| -13.5| 3.80|   -1813|    3.06|  0.000|
|2          |2_8861444_A_E    |   8861444|A      |E   |      0|E       |A       | 0.127| -14.1| 3.76|   -1813|    3.18|  0.000|
|2          |2_9196139_H_E    |   9196139|H      |E   |      0|E       |H       | 0.128| -12.7| 3.81|   -1814|    3.04|  0.001|
|5          |5_117508066_B_E  | 117508066|B      |E   |      0|E       |B       | 0.153|  13.2| 3.58|   -1813|    2.80|  0.000|
|5          |5_123629774_B_E  | 123629774|B      |E   |      0|E       |B       | 0.099|  16.4| 4.36|   -1812|    3.00|  0.000|
|5          |5_39679307_E_A   |  39679307|E      |A   |      0|A       |E       | 0.090| -15.7| 4.39|   -1813|    2.83|  0.000|
|5          |5_40272493_E_A   |  40272493|E      |A   |      0|A       |E       | 0.095| -15.6| 4.30|   -1813|    2.79|  0.000|
|5          |5_54858991_E_A   |  54858991|E      |A   |      0|A       |E       | 0.102| -14.6| 4.24|   -1814|    2.87|  0.001|
|5          |5_55215630_E_A   |  55215630|E      |A   |      0|A       |E       | 0.101| -14.7| 4.26|   -1814|    2.88|  0.001|
|6          |6_107008934_D_A  | 107008934|D      |A   |      0|A       |D       | 0.077| -16.9| 4.56|   -1813|    3.23|  0.000|
|8          |8_3000000_H_C    |   3000000|H      |C   |      0|C       |H       | 0.077| -15.1| 4.35|   -1814|    3.00|  0.001|
|8          |8_3410751_H_C    |   3410751|H      |C   |      0|C       |H       | 0.077| -15.2| 4.36|   -1813|    3.01|  0.001|
|10         |10_100114175_D_D | 100114175|D      |D   |      0|D       |D       | 0.171| -11.0| 3.25|   -1814|    2.30|  0.001|
|10         |10_57753923_G_C  |  57753923|G      |C   |      0|C       |G       | 0.101|  15.3| 4.38|   -1813|    2.66|  0.001|
|10         |10_58072086_G_C  |  58072086|G      |C   |      0|C       |G       | 0.101|  15.1| 4.39|   -1814|    2.66|  0.001|
|10         |10_60353971_F_C  |  60353971|F      |C   |      0|C       |F       | 0.084|  16.4| 4.78|   -1814|    2.71|  0.001|
|10         |10_95418856_F_E  |  95418856|F      |E   |      0|E       |F       | 0.120|  13.1| 3.86|   -1814|    3.28|  0.001|
|10         |10_96637576_F_E  |  96637576|F      |E   |      0|E       |F       | 0.116|  13.1| 3.92|   -1814|    3.28|  0.001|
|10         |10_99490238_H_D  |  99490238|H      |D   |      0|D       |H       | 0.163| -12.6| 3.26|   -1812|    2.17|  0.000|
|12         |12_108748974_D_C | 108748974|D      |C   |      0|C       |D       | 0.122|  14.9| 3.78|   -1812|    2.92|  0.000|
|12         |12_110011843_A_C | 110011843|A      |C   |      0|C       |A       | 0.121|  12.9| 3.86|   -1814|    2.96|  0.001|
|13         |13_26871268_G_C  |  26871268|G      |C   |      0|C       |G       | 0.147|  12.4| 3.67|   -1814|    2.57|  0.001|
|13         |13_27924769_G_C  |  27924769|G      |C   |      0|C       |G       | 0.147|  13.3| 3.64|   -1813|    2.63|  0.000|
|13         |13_28620475_G_C  |  28620475|G      |C   |      0|C       |G       | 0.153|  14.0| 3.56|   -1812|    2.56|  0.000|
|13         |13_29976363_E_C  |  29976363|E      |C   |      0|C       |E       | 0.161|  15.5| 3.49|   -1810|    2.57|  0.000|
|13         |13_30180778_E_C  |  30180778|E      |C   |      0|C       |E       | 0.161|  15.5| 3.49|   -1810|    2.57|  0.000|
|13         |13_31072103_E_C  |  31072103|E      |C   |      0|C       |E       | 0.159|  13.9| 3.52|   -1812|    2.56|  0.000|
|13         |13_32418268_A_C  |  32418268|A      |C   |      0|C       |A       | 0.155|  13.0| 3.54|   -1813|    2.49|  0.000|
|16         |16_63647011_B_B  |  63647011|B      |B   |      0|B       |B       | 0.156| -12.2| 3.60|   -1814|    3.60|  0.001|
|19         |19_25487123_G_G  |  25487123|G      |G   |      0|G       |G       | 0.129|  14.2| 3.98|   -1813|    3.05|  0.000|

``` r
manhattan_plot(x = lmm.hfd.results, pval.colname = "p_wald", chr.colname = "chromosome", pos.colname = "position", plot.title = "DO Mice on a HFHS Diet", y.label = "LOD Score") -> hfd.manhattan

hfd.manhattan
```

![](figures/chol-snp-analysis-hfd-2.png)<!-- -->

#### Chromosome 13 QTL Associated with Cholesterol on HFHS Diets


``` r
snp.pos <- 29976363
library(forcats)
peak13 <- dplyr::filter(lmm.hfd.results,
                  chromosome==13,
                  position>snp.pos-10000000,
                  position<snp.pos+10000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H")) %>%
  mutate(LOD=-log10(p_wald)) %>%
  mutate(LOD.drop = LOD-max(LOD))
library(ggplot2)
ggplot(data=peak13,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size on HFHS (mg/dL)",
       x="Position on Chromosome 13") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=16) +
  theme(legend.position=c(0.8,0.2),
        legend.text=element_text(size=8))
```

![](figures/peak-13-hfd-1.png)<!-- -->

``` r
#credible region 1.5 LOD less than peak 28011631-32375634
#From GenomMUSTer 3755 variants differ btween 129 and bl6
#upstream variants in Sox4 and 
#missense variant in rs29850511	chr13:29326051	T	Cdkal1	ENSMUST00000006353.13	Transcript	missense_variant	1879	1723	575	L/M	Ttg/Atg	rs29850511	EXON=16/16
#GWAS associated Cdkal1 with HDL and TC, seems to be a negative regulator of bile acid production, is atheroprotective in Apoe null mice
```

#### Chromosome 18 QTL Associated with Cholesterol on NCD Diets


``` r
snp.pos <- 46410922
library(forcats)
peak13 <- dplyr::filter(lmm.ncd.results,
                  chromosome==18,
                  position>snp.pos-10000000,
                  position<snp.pos+10000000) %>%
  mutate(alt=fct_recode(as.factor(alt),
                        "C57BL/6J"="A",
                        "NZO"="B",
                        "Sv129"="C",
                        "PWK"="D",
                        "A/J"="E",
                        "NOD"="F",
                        "CAST"="G",
                        "WSB"="H")) %>%
  mutate(LOD=-log10(p_wald)) %>%
  mutate(LOD.drop = LOD-max(LOD))
library(ggplot2)
ggplot(data=peak13,
       aes(x=position,
       y=beta,
       col=alt,
       group=alt)) +
  geom_line() +
  geom_hline(yintercept=0,lty=2) +
  labs(title="Strain Specific Effects",
       y="Cholesterol Effect Size on NCD (mg/dL)",
       x="Position on Chromosome 18") +
  scale_color_discrete(name="") +
  guides(col=guide_legend(ncol=2)) +
  theme_classic(base_size=16) +
  theme(legend.position=c(0.8,0.2),
        legend.text=element_text(size=8))
```

![](figures/peak-18-ncd-1.png)<!-- -->

``` r
#credible region 1.5 LOD less than peak 28011631-32375634
#From GenomMUSTer 3755 variants differ btween 129 and bl6
#upstream variants in Sox4 and 
#missense variant in rs29850511	chr13:29326051	T	Cdkal1	ENSMUST00000006353.13	Transcript	missense_variant	1879	1723	575	L/M	Ttg/Atg	rs29850511	EXON=16/16
#GWAS associated Cdkal1 with HDL and TC, seems to be a negative regulator of bile acid production, is atheroprotective in Apoe null mice
```

### Comparason of NCD and HFHS GWAS


``` r
library(cowplot)

# plots are drawn without alignment
plot_grid(ncd.manhattan, hfd.manhattan, align="v",ncol=1)
```

![](figures/ncd-hfd-comparason-1.png)<!-- -->

# Summary of Interesting SNPs


``` r
snps.of.interest <- bind_rows(additive.snp.summary,hfhs.snp.summary,ncd.snp.summary)

write_csv(snps.of.interest, file="SNPs_of_interest_diet.csv")
snps.of.interest %>% distinct(rs) %>% pull(rs) %>% write(file="SNPs_of_interest_diet.txt") #just the rsids
```

# Session Information


``` r
sessionInfo()
```

```
## R version 4.4.0 (2024-04-24)
## Platform: x86_64-pc-linux-gnu
## Running under: Red Hat Enterprise Linux 8.8 (Ootpa)
## 
## Matrix products: default
## BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRblas.so 
## LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: America/Detroit
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] cowplot_1.1.3 forcats_1.0.0 ggmanh_1.8.0  ggplot2_3.5.1 qqman_0.1.9  
##  [6] broom_1.0.6   dplyr_1.1.4   tidyr_1.3.1   readr_2.1.5   knitr_1.48   
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.9         utf8_1.2.4         generics_0.1.3     hms_1.1.3         
##  [5] digest_0.6.36      magrittr_2.0.3     evaluate_0.24.0    grid_4.4.0        
##  [9] RColorBrewer_1.1-3 calibrate_1.7.7    fastmap_1.2.0      jsonlite_1.8.8    
## [13] backports_1.5.0    purrr_1.0.2        fansi_1.0.6        scales_1.3.0      
## [17] jquerylib_0.1.4    cli_3.6.3          rlang_1.1.4        crayon_1.5.3      
## [21] bit64_4.0.5        munsell_0.5.1      withr_3.0.0        cachem_1.1.0      
## [25] yaml_2.3.9         tools_4.4.0        parallel_4.4.0     tzdb_0.4.0        
## [29] colorspace_2.1-0   vctrs_0.6.5        R6_2.5.1           lifecycle_1.0.4   
## [33] bit_4.0.5          vroom_1.6.5        MASS_7.3-60.2      pkgconfig_2.0.3   
## [37] pillar_1.9.0       bslib_0.7.0        gtable_0.3.5       glue_1.7.0        
## [41] highr_0.11         xfun_0.45          tibble_3.2.1       tidyselect_1.2.1  
## [45] rstudioapi_0.16.0  farver_2.1.2       htmltools_0.5.8.1  rmarkdown_2.27    
## [49] labeling_0.4.3     compiler_4.4.0
```

