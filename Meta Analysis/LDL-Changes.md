---
title: "Systematic Review of LDL Changes on a Ketogenic Diet"
author: "Cody Cousineau, Kai-Lin Jen, and Dave Bridges"
date: "June 18, 2021"
output:
  html_document:
    highlight: tango
    keep_md: yes
    number_sections: yes
    toc: yes
  pdf_document:
    toc: yes
---


# Purpose

To evaluate how LDL changes in ketogenic diet studies

# Experimental Details

Evaluated studies where ketogenic diets (<25g/day of CHO) are used and weight and LDL-C are reported as an outcome

# Raw Data

Reviewed data from the Choi *et al* meta-analysis (http://dx.doi.org/10.3390/nu12072005), pulling in data on baseline weight, weight changes, LDL, LDL changes and standard deviations. A systematic literature search of PubMed was then performed to identify other randomized controlled trials (RCTs) and single-arm interventions of patients that evaluated the effects of a ketogenic diet on weight and lipid profile as primary endpoints. All studies using a KD diet that met our inclusion criteria where intake of carbohydrate was less than 25 grams per day were included. This search was most recently updated on Fri Mar 15 09:58:07 2024.

We used a value 130mg/dL of LDL-C at baseline to stratify individuals as being hypercholesterolemic or not.


Correlations in this meta-analysis before and after the administration of the ketogenic diet were analyzed with linear models and results given using pearsons correlation coefficient, statistical significance was defined as below 0.05. 


For all outcomes, we tested sex as a modifier and as a covariate. For outcomes where sex was found to be a significant modifier, these results are reported. 


```r
filename <- 'VLCF Meta-Analysis.csv'
#filename <- 'LDL Study Summary.xlsx' #make this a separate line, you can use any variable you want
google.sheet.link <- 'https://docs.google.com/spreadsheets/d/e/2PACX-1vRkVL0hHHPMJ_fI8EVc_RiFJvcnTHmBySUd3bvikUJjancG8DCsaG5k0eIFUBfov4drTr8MMNV8GLnv/pub?gid=0&single=true&output=csv'
download.file(google.sheet.link,destfile=filename) #can commnt this out if you dont want to udpat the file

#this loads whatever the file is into a dataframe called exp.data if it exists
#make a new variable called break.groups, make it so c is the break.groups so that is the code

exp.data <- read_csv(filename)
eval.data <- 
  exp.data %>% 
  filter(Use=='x') %>%
  mutate(across(`Baseline Weight`:n,.fns=as.numeric)) %>% #force all fields to be numeric
  mutate(Pct.Wt.Change = `Weight Change`/`Baseline Weight`*100)%>%
  mutate(PCT.BMI.Change = `BMI Change` / `Baseline BMI`*100)%>%
  mutate(Sex.Group = cut(`Percent Male`, breaks = c(0,.1,.9,1), include.lowest = TRUE, labels = c("Mostly Female", "Mixed", "Mostly Male")))
```

These data can be found in **/Users/davebrid/Documents/GitHub/PrecisionNutrition/Meta Analysis** in a file named **VLCF Meta-Analysis.csv**.  This script was most recently updated on **Fri Mar 15 09:58:09 2024**.

# Meta-Analysis 


```r
meta.data <-
  eval.data %>%
  filter(!is.na(`LDL Endpoint SD`)) %>%
  mutate(Pooled.SD=sqrt(`Baseline SD`^2+`LDL Endpoint SD`^2)) %>%
  mutate(SMD=`Change in LDL-C`) %>%
  mutate(SMD.Wt = `Weight Change`) %>%
  mutate(BMD=`BMI Change`)

library(meta)
ldl.c.meta <- metagen(TE = SMD,
                 seTE = Pooled.SD,
                 n.e=n,
                 n.c=n,
                 studlab = Study,
                 data = meta.data,
                 sm = "SMD",
                 comb.fixed = TRUE,
                 comb.random = FALSE,
                 #method.tau = "REML",
                 hakn = TRUE,
                 title = "LDL-C Changes in Ketogenic Diet Studies")


forest(ldl.c.meta, 
            sortvar = SMD,
            #predict = F, 
            print.tau2 = TRUE,
            print.I2 = TRUE,
            #leftcols=c('Study'),
            #rightlabs=c('Change','95% CI','Weight'),
           layout = "RevMan5",
          fontsize=8,
          leftcols = c("Study",'effect',"seTE",'ci','w.fixed'),
          leftlabs = c("Study","Change","SE","95% CI","Weight"),
          range = F,
          digits=1,
          digits.se=1)
```

![](figures/meta-analysis-1.png)<!-- -->

```r
bmi.meta <- metagen(TE = BMD,
                 seTE = Pooled.SD,
                 n.e=n,
                 n.c=n,
                 studlab = Study,
                 data = meta.data,
                 sm = "BMD",
                 comb.fixed = TRUE,
                 comb.random = FALSE,
                 #method.tau = "REML",
                 hakn = TRUE,
                 title = "Weight Changes in Ketogenic Diet Studies")

forest(bmi.meta, 
            sortvar = BMD,
            #predict = F, 
            print.tau2 = TRUE,
            print.I2 = TRUE,
            #leftcols=c('Study'),
            #rightlabs=c('Change','95% CI','Weight'),
           layout = "RevMan5",
          fontsize=8,
          leftcols = c("Study",'effect',"seTE",'ci','w.fixed'),
          leftlabs = c("Study","Change","SE","95% CI","Weight"),
          range = F,
          digits=1,
          digits.se=1)
```

![](figures/meta-analysis-2.png)<!-- -->

We evaluated 2 studies for this meta-analysis. Using the meta-analysis method, we found fasting blood LDL-C levels were increased -2.1 mg/dL (95% CI: -59.378 to 55.178) after the ketogenic diet intervention compared to pre-intervention levels, with a significant p-value of 0.943. Across these studies, the I<sup>2</sup> is 0, the p-value for Q is 1. This is a highly consistent I^2. 

# Average Change in LDL-C


```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=reorder(Study, -`Change in LDL-C`))) +
  geom_bar(stat='identity') +
  labs(x="Study")
```

![](figures/ldl-change-1.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=reorder(Study, -`Change in LDL-C`))) +
  geom_bar(stat='identity', aes(fill=`Percent Male`)) +
  labs(x="Study")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](figures/ldl-change-2.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=reorder(Study, -`Change in LDL-C`))) +
  geom_bar(stat='identity', aes(fill=`Sex.Group`)) +
  labs(x="Study") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](figures/ldl-change-3.png)<!-- -->

## Relative to Weight


```r
library(ggrepel)
library(ggplot2)
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline Weight`)) +
  geom_point() +
  geom_smooth(span=1) +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline Weight")
```

![](figures/ldl-change-vs-weight-1.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline BMI`)) +
  geom_point() +
  geom_smooth(span=1) +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline BMI")
```

![](figures/ldl-change-vs-weight-2.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline Weight`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline Weight")
```

![](figures/ldl-change-vs-weight-3.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline BMI`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline BMI")
```

![](figures/ldl-change-vs-weight-4.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline Weight`,
             col=`Percent Male`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline Weight")
```

![](figures/ldl-change-vs-weight-5.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline Weight`,
             col=`Sex.Group`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline Weight")
```

![](figures/ldl-change-vs-weight-6.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline BMI`,
             col=`Sex.Group`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline Weight")
```

![](figures/ldl-change-vs-weight-7.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline BMI`,
             col=`Normal weight`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline Weight")
```

![](figures/ldl-change-vs-weight-8.png)<!-- -->

```r
library(broom)
bind_rows(`Baseline`=shapiro.test(eval.data$`Baseline Weight`)$p.value,
          `Change`=shapiro.test(eval.data$`Change in LDL-C`)$p.value) %>%
  kable(caption="Shapiro Tests for Correlates")
```



Table: Shapiro Tests for Correlates

| Baseline| Change|
|--------:|------:|
|    0.233|      0|

```r
with(eval.data, cor.test(`Change in LDL-C`,`Baseline Weight`, method="spearman")) %>% tidy %>% kable(caption="Correlation between change in LDL-C and baseline weight")
```



Table: Correlation between change in LDL-C and baseline weight

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|   -0.765|      2347|       0|Spearman's rank correlation rho |two.sided   |

```r
lm(`Change in LDL-C`~`Baseline Weight`+`Percent Male`, data=eval.data) %>% tidy %>% kable(caption="Linear model between change in LDL-C and baseline weight, including gender")
```



Table: Linear model between change in LDL-C and baseline weight, including gender

|term              | estimate| std.error| statistic| p.value|
|:-----------------|--------:|---------:|---------:|-------:|
|(Intercept)       |    70.52|    15.729|     4.483|   0.000|
|`Baseline Weight` |    -0.64|     0.196|    -3.263|   0.005|
|`Percent Male`    |    -4.11|    12.169|    -0.337|   0.740|

```r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight`*`Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term                        | df| sumsq| meansq| statistic| p.value|
|:---------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight`           |  1|  3436|   3436|    16.075|   0.001|
|Sex.Group                   |  2|   224|    112|     0.524|   0.603|
|`Baseline Weight`:Sex.Group |  2|   309|    154|     0.722|   0.503|
|Residuals                   | 14|  2993|    214|        NA|      NA|

```r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  3436|   3436|    16.653|   0.001|
|Sex.Group         |  2|   224|    112|     0.543|   0.592|
|Residuals         | 16|  3302|    206|        NA|      NA|

```r
ldl.weight.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline Weight`, data = eval.data)
ldl.weight.baseline.lm %>% tidy %>% kable
```



|term              | estimate| std.error| statistic| p.value|
|:-----------------|--------:|---------:|---------:|-------:|
|(Intercept)       |   72.153|    14.589|      4.95|   0.000|
|`Baseline Weight` |   -0.676|     0.161|     -4.19|   0.001|

```r
ldl.weightsex.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |    67.53|    15.949|     4.235|   0.001|
|`Baseline Weight`    |    -0.56|     0.213|    -2.636|   0.018|
|Sex.GroupMixed       |    -7.12|     9.166|    -0.776|   0.449|
|Sex.GroupMostly Male |   -12.96|    12.875|    -1.006|   0.329|

```r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  3436|   3436|    16.653|   0.001|
|Sex.Group         |  2|   224|    112|     0.543|   0.592|
|Residuals         | 16|  3302|    206|        NA|      NA|

```r
ldl.bmi.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, data = eval.data)
ldl.bmi.baseline.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    57.41|    14.110|      4.07|   0.002|
|`Baseline BMI` |    -1.52|     0.444|     -3.43|   0.006|

```r
ldl.bmi.baseline.sex.aov <- aov(`Change in LDL-C` ~ `Baseline BMI`*`Sex.Group`, data = eval.data)
ldl.bmi.baseline.sex.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI`           |  1|   920|  920.2|     29.39|   0.001|
|Sex.Group                |  2|   393|  196.6|      6.28|   0.023|
|`Baseline BMI`:Sex.Group |  1|   218|  217.9|      6.96|   0.030|
|Residuals                |  8|   251|   31.3|        NA|      NA|

```r
ldl.bmi.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline BMI` + `Sex.Group`, data = eval.data)
ldl.bmi.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI` |  1|   920|    920|     17.68|   0.002|
|Sex.Group      |  2|   393|    197|      3.78|   0.064|
|Residuals      |  9|   468|     52|        NA|      NA|

```r
ldl.bmi.baseline.mixed.lm<- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mixed"))
ldl.bmi.baseline.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mixed groups")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mixed groups

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   20.410|    14.085|      1.45|   0.191|
|`Baseline BMI` |   -0.478|     0.415|     -1.15|   0.288|

```r
ldl.bmi.baseline.female.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mostly Female"))
ldl.bmi.baseline.female.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mostly Female groups")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mostly Female groups

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    95.07|    21.041|      4.52|   0.139|
|`Baseline BMI` |    -2.61|     0.781|     -3.34|   0.185|

```r
ldl.bmi.baseline.male.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mostly Male"))
ldl.bmi.baseline.male.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mostly Male")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mostly Male

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |     10.7|       NaN|       NaN|     NaN|
|`Baseline BMI` |       NA|        NA|        NA|      NA|

```r
ldl.bmi.baseline.sex.aov <- aov(`Change in LDL-C` ~ `Baseline BMI`*`Normal weight`, data = eval.data)
ldl.bmi.baseline.sex.aov %>% tidy %>% kable
```



|term                           | df| sumsq| meansq| statistic| p.value|
|:------------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI`                 |  1| 920.2|  920.2|     9.877|   0.020|
|`Normal weight`                |  4| 214.4|   53.6|     0.575|   0.692|
|`Baseline BMI`:`Normal weight` |  1|  88.2|   88.2|     0.947|   0.368|
|Residuals                      |  6| 559.0|   93.2|        NA|      NA|

```r
ldl.bmi.baseline.female.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, `Normal weight` == "Yes"))
ldl.bmi.baseline.female.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Normal weight studies")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Normal weight studies

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    109.1|    137.49|     0.793|   0.511|
|`Baseline BMI` |     -3.6|      5.75|    -0.626|   0.595|

```r
ldl.bmi.baseline.male.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, `Normal weight` == "No"))
ldl.bmi.baseline.male.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Overweight studies")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Overweight studies

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |  -14.311|     35.50|    -0.403|   0.707|
|`Baseline BMI` |    0.453|      1.02|     0.446|   0.679|

Lower baseline BMI was associated with an increased change in LDL-C after consumption of a ketogenic diet (r<sup>2</sup> = 0.516, p-value = 0.006). The association with increased LDL-C was consistent with baseline weight, where a lower baseline weight was associated with an increased change in LDL-C (r<sup>2</sup> = 0.494, p-value = 0.001). 

## Relative to Weight Loss


```r
library(ggplot2)

eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=Pct.Wt.Change)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Weight Change (%)")
```

![](figures/ldl-change-vs-weight-change-1.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`BMI Change`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="BMI Change")
```

![](figures/ldl-change-vs-weight-change-2.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`BMI Change`,
             col = `Sex.Group`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="BMI Change")
```

![](figures/ldl-change-vs-weight-change-3.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`PCT.BMI.Change`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="% BMI Change")
```

![](figures/ldl-change-vs-weight-change-4.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`PCT.BMI.Change`,
             col = `Sex.Group`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="% BMI Change")
```

![](figures/ldl-change-vs-weight-change-5.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=Pct.Wt.Change,
             col=`Percent Male`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Weight Change (%)")
```

![](figures/ldl-change-vs-weight-change-6.png)<!-- -->

```r
library(broom)
bind_rows(`Weight Change`=shapiro.test(eval.data$Pct.Wt.Change)$p.value,
          `LDL Change`=shapiro.test(eval.data$`Change in LDL-C`)$p.value) %>%
  kable(caption="Shapiro Tests for Correlates")
```



Table: Shapiro Tests for Correlates

| Weight Change| LDL Change|
|-------------:|----------:|
|         0.266|          0|

```r
with(eval.data, cor.test(`Change in LDL-C`,Pct.Wt.Change, method="spearman")) %>% tidy %>% kable(caption="Correlation between change in LDL-C and weight change")
```



Table: Correlation between change in LDL-C and weight change

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|    0.476|       598|   0.039|Spearman's rank correlation rho |two.sided   |

```r
ldl.weightsexinter.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change`*`Sex.Group`, data = eval.data)
ldl.weightsexinter.change.aov %>% tidy %>% kable
```



|term                    | df| sumsq| meansq| statistic| p.value|
|:-----------------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change           |  1|   672|  671.9|     4.021|   0.066|
|Sex.Group               |  2|   452|  226.0|     1.352|   0.293|
|Pct.Wt.Change:Sex.Group |  2|   115|   57.7|     0.345|   0.714|
|Residuals               | 13|  2173|  167.1|        NA|      NA|

```r
ldl.weightsex.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change` + `Sex.Group`, data = eval.data)
ldl.weightsex.change.aov %>% tidy %>% kable
```



|term          | df| sumsq| meansq| statistic| p.value|
|:-------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change |  1|   672|    672|      4.41|   0.053|
|Sex.Group     |  2|   452|    226|      1.48|   0.259|
|Residuals     | 15|  2288|    153|        NA|      NA|

```r
ldl.weight.change.lm <- lm(`Change in LDL-C` ~ `Pct.Wt.Change`, data = eval.data)
ldl.weight.change.lm %>% tidy %>% kable
```



|term          | estimate| std.error| statistic| p.value|
|:-------------|--------:|---------:|---------:|-------:|
|(Intercept)   |    18.80|      5.45|      3.45|   0.003|
|Pct.Wt.Change |     1.51|      0.74|      2.04|   0.057|

```r
ldl.bmisexinter.change.aov <- aov(`Change in LDL-C` ~ `BMI Change`*`Sex.Group`, data = eval.data)
ldl.bmisexinter.change.aov %>% tidy %>% kable
```



|term                   | df| sumsq| meansq| statistic| p.value|
|:----------------------|--:|-----:|------:|---------:|-------:|
|`BMI Change`           |  1|   674|  674.1|     15.61|   0.011|
|Sex.Group              |  2|   496|  248.0|      5.74|   0.051|
|`BMI Change`:Sex.Group |  1|   380|  380.1|      8.80|   0.031|
|Residuals              |  5|   216|   43.2|        NA|      NA|

```r
ldl.bmisex.change.aov <- aov(`Change in LDL-C` ~ `BMI Change` + `Sex.Group`, data = eval.data)
ldl.bmisex.change.aov %>% tidy %>% kable
```



|term         | df| sumsq| meansq| statistic| p.value|
|:------------|--:|-----:|------:|---------:|-------:|
|`BMI Change` |  1|   674|  674.1|      6.79|   0.040|
|Sex.Group    |  2|   496|  248.0|      2.50|   0.163|
|Residuals    |  6|   596|   99.3|        NA|      NA|

```r
ldl.bmi.change.lm <- lm(`Change in LDL-C` ~ `BMI Change`, data = eval.data)
ldl.bmi.change.lm %>% tidy %>% kable
```



|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |    23.16|      6.97|      3.32|   0.011|
|`BMI Change` |     5.57|      2.51|      2.22|   0.057|

```r
ldl.bmi.change.mixed.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mixed"))
ldl.bmi.change.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mixed groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mixed groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |    1.168|      6.56|     0.178|   0.867|
|`BMI Change` |   -0.261|      1.90|    -0.138|   0.897|

```r
ldl.bmi.change.male.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Male"))
ldl.bmi.change.male.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Male groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mostly Male groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |     10.7|       NaN|       NaN|     NaN|
|`BMI Change` |       NA|        NA|        NA|      NA|

```r
ldl.bmi.change.female.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Female"))
ldl.bmi.change.female.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Female groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mostly Female groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |     46.8|     10.40|      4.50|   0.139|
|`BMI Change` |     19.1|      8.25|      2.31|   0.260|

```r
ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
ldl.pctbmi.change.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    22.70|     7.868|      2.89|    0.02|
|PCT.BMI.Change |     1.75|     0.939|      1.86|    0.10|

```r
ldl.pctbmisexinter.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change`*`Sex.Group`, data = eval.data)
ldl.pctbmisexinter.change.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|PCT.BMI.Change           |  1|   533|  533.0|      9.40|   0.028|
|Sex.Group                |  2|   615|  307.5|      5.42|   0.056|
|PCT.BMI.Change:Sex.Group |  1|   335|  334.8|      5.91|   0.059|
|Residuals                |  5|   283|   56.7|        NA|      NA|

```r
ldl.pctbmisex.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change` + `Sex.Group`, data = eval.data)
ldl.pctbmisex.change.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|PCT.BMI.Change |  1|   533|    533|      5.17|   0.063|
|Sex.Group      |  2|   615|    307|      2.98|   0.126|
|Residuals      |  6|   618|    103|        NA|      NA|

```r
ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
ldl.pctbmi.change.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    22.70|     7.868|      2.89|    0.02|
|PCT.BMI.Change |     1.75|     0.939|      1.86|    0.10|

```r
ldl.pctbmi.change.mixed.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mixed"))
ldl.pctbmi.change.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mixed groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mixed groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |    1.168|      6.56|     0.178|   0.867|
|`BMI Change` |   -0.261|      1.90|    -0.138|   0.897|

```r
ldl.pctbmi.change.male.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Male"))
ldl.pctbmi.change.male.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Male groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mostly Male groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |     10.7|       NaN|       NaN|     NaN|
|`BMI Change` |       NA|        NA|        NA|      NA|

```r
ldl.pctbmi.change.female.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Female"))
ldl.pctbmi.change.female.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Female groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mostly Female groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |     46.8|     10.40|      4.50|   0.139|
|`BMI Change` |     19.1|      8.25|      2.31|   0.260|
Greater BMI decreases over the study period were associated with a smaller increase in LDL-C after consumption of a ketogenic diet, though this did not reach significance (r<sup>2</sup> = 0.382, p-value = 0.057). The association with the change in LDL-C and decrease in BMI was consistent with weight, with change in weight on LDL-C reaching significance (r<sup>2</sup> = 0.197, p-value = 0.057), where greater decreases in weight were associated with lower increases in LDL-C after consumption of a ketogenic diet. Looking at percent BMI change to account for baseline BMI, greater percent change decreases were associated with a lower increase in LDL-C on a ketogenic diet, though this was not significant (r<sup>2</sup> = 0.302, p-value = 0.1).


## Relative to Baseline LDL-C


```r
library(ggplot2)
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline LDL`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline LDL-C")
```

![](figures/ldl-change-vs-baseline-1.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline LDL`,
             col=`Percent Male`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline LDL-C")
```

![](figures/ldl-change-vs-baseline-2.png)<!-- -->

```r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline LDL`,
             col=`Sex.Group`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline LDL-C")
```

![](figures/ldl-change-vs-baseline-3.png)<!-- -->

```r
library(broom)
bind_rows(`Baseline`=shapiro.test(eval.data$`Baseline LDL`)$p.value,
          `Change`=shapiro.test(eval.data$`Change in LDL-C`)$p.value) %>%
  kable(caption="Shapiro Tests for Correlates")
```



Table: Shapiro Tests for Correlates

| Baseline| Change|
|--------:|------:|
|    0.258|      0|

```r
with(eval.data, cor.test(`Change in LDL-C`,`Baseline LDL`, method="spearman")) %>% tidy %>% kable(caption="Correlation between baseline and change in LDL-C")
```



Table: Correlation between baseline and change in LDL-C

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|   -0.299|      2300|   0.177|Spearman's rank correlation rho |two.sided   |

```r
ldl.baselinesexinter.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` * `Sex.Group`, data = eval.data)
ldl.baselinesexinter.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline LDL`           |  1|   794|    794|      3.94|   0.064|
|Sex.Group                |  2|  2304|   1152|      5.72|   0.013|
|`Baseline LDL`:Sex.Group |  2|   894|    447|      2.22|   0.141|
|Residuals                | 16|  3222|    201|        NA|      NA|

```r
ldl.baselinesex.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baselinesex.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline LDL` |  1|   794|    794|      3.47|   0.079|
|Sex.Group      |  2|  2304|   1152|      5.04|   0.018|
|Residuals      | 18|  4117|    229|        NA|      NA|

```r
ldl.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline LDL`, data = eval.data)
ldl.baseline.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   43.404|    20.678|      2.10|   0.049|
|`Baseline LDL` |   -0.289|     0.184|     -1.57|   0.132|

```r
ldl.baselinesex.lm <- lm(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baselinesex.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |   58.435|    18.316|      3.19|   0.005|
|`Baseline LDL`       |   -0.277|     0.159|     -1.74|   0.099|
|Sex.GroupMixed       |  -21.066|     7.529|     -2.80|   0.012|
|Sex.GroupMostly Male |  -28.078|    10.767|     -2.61|   0.018|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline LDL` |  1|   794|    794|      3.47|   0.079|
|Sex.Group      |  2|  2304|   1152|      5.04|   0.018|
|Residuals      | 18|  4117|    229|        NA|      NA|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term      | df| sumsq| meansq| statistic| p.value|
|:---------|--:|-----:|------:|---------:|-------:|
|Sex.Group |  2|  2406|   1203|      4.75|   0.021|
|Residuals | 19|  4809|    253|        NA|      NA|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Percent Male`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Percent Male` |  1|  1266|   1266|      4.26|   0.052|
|Residuals      | 20|  5949|    297|        NA|      NA|

Among individuals, baseline LDL-C was not positively correlated with change in LDL-C after consumption of a ketogenic diet and the relationship was not significant (r<sup>2</sup> = 0.11, p-value = 0.132). 

# Session Information


```r
sessionInfo()
```

```
## R version 4.2.2 (2022-10-31)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] broom_1.0.5   ggrepel_0.9.5 meta_7.0-0    metadat_1.2-0 readr_2.1.5  
## [6] ggplot2_3.5.0 dplyr_1.1.4   tidyr_1.3.1   knitr_1.45   
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.12         lattice_0.22-5      digest_0.6.35      
##  [4] utf8_1.2.4          R6_2.5.1            backports_1.4.1    
##  [7] evaluate_0.23       highr_0.10          pillar_1.9.0       
## [10] rlang_1.1.3         rstudioapi_0.15.0   minqa_1.2.6        
## [13] jquerylib_0.1.4     nloptr_2.0.3        Matrix_1.5-4.1     
## [16] rmarkdown_2.26      mathjaxr_1.6-0      labeling_0.4.3     
## [19] splines_4.2.2       lme4_1.1-35.1       stringr_1.5.1      
## [22] bit_4.0.5           munsell_0.5.0       compiler_4.2.2     
## [25] numDeriv_2016.8-1.1 xfun_0.42           pkgconfig_2.0.3    
## [28] mgcv_1.9-1          htmltools_0.5.7     tidyselect_1.2.1   
## [31] tibble_3.2.1        fansi_1.0.6         crayon_1.5.2       
## [34] tzdb_0.4.0          withr_3.0.0         MASS_7.3-60.0.1    
## [37] grid_4.2.2          nlme_3.1-164        jsonlite_1.8.8     
## [40] gtable_0.3.4        lifecycle_1.0.4     magrittr_2.0.3     
## [43] metafor_4.4-0       scales_1.3.0        cli_3.6.2          
## [46] stringi_1.8.3       vroom_1.6.5         cachem_1.0.8       
## [49] farver_2.1.1        xml2_1.3.6          bslib_0.6.1        
## [52] generics_0.1.3      vctrs_0.6.5         boot_1.3-30        
## [55] tools_4.2.2         bit64_4.0.5         CompQuadForm_1.4.3 
## [58] glue_1.7.0          purrr_1.0.2         hms_1.1.3          
## [61] parallel_4.2.2      fastmap_1.1.1       yaml_2.3.8         
## [64] colorspace_2.1-0    sass_0.4.8
```
