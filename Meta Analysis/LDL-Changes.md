---
title: "Systematic Review of LDL Changes on a Ketogenic Diet"
author: "Cody Cousineau and Dave Bridges"
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

Reviewed data from the Choi *et al* meta-analysis (http://dx.doi.org/10.3390/nu12072005), pulling in data on baseline weight, weight changes, LDL, LDL changes and standard deviations. A systematic literature search of PubMed was then performed to identify other randomized controlled trials (RCTs) and single-arm interventions of patients that evaluated the effects of a ketogenic diet on weight and lipid profile as primary endpoints. All studies using a KD diet that met our inclusion criteria where intake of carbohydrate was less than 25 grams per day were included. This search was most recently updated on Thu Jan 11 13:36:08 2024.

We used a value 130mg/dL of LDL-C at baseline to stratify individuals as being hypercholesterolemic or not.


Correlations in this meta-analysis before and after the administration of the ketogenic diet were analyzed with linear models and results given using pearsons correlation coefficient, statistical significance was defined as below 0.05. 


For all outcomes, we tested sex as a modifier and as a covariate. For outcomes where sex was found to be a significant modifier, these results are reported. 


```r
library(readxl) #loads the readr package
filename <- 'LDL Study Summary.xlsx' #make this a separate line, you can use any variable you want

#this loads whatever the file is into a dataframe called exp.data if it exists
#make a new variable called break.groups, make it so c is the break.groups so that is the code

exp.data <- read_excel(filename)
eval.data <- exp.data %>% 
  filter(Use=='x') %>%
  mutate(Pct.Wt.Change = `Weight Change`/`Baseline Weight`*100)%>%
  mutate(PCT.BMI.Change = `BMI Change` / `Baseline BMI`*100)%>%
  mutate(Sex.Group = cut(`Percent Male`, breaks = c(0,.1,.9,1), include.lowest = TRUE, labels = c("Mostly Female", "Mixed", "Mostly Male")))
```

These data can be found in **/Users/davebrid/Documents/GitHub/PrecisionNutrition/Meta Analysis** in a file named **LDL Study Summary.xlsx**.  This script was most recently updated on **Thu Jan 11 13:36:08 2024**.

# Analysis

# Variance in LDL-C


```r
library(ggplot2)

eval.data <-
  eval.data %>%
  mutate(LDL.endpoint = `Baseline LDL`+`Change in LDL-C`)

eval.data %>%
  filter(!(is.na(`Endpoint SD`))) %>%
  ggplot(aes(y=LDL.endpoint,
             ymin=LDL.endpoint-`Endpoint SD`,
             ymax=LDL.endpoint+`Endpoint SD`,
             x=reorder(Study,-`LDL.endpoint` ))) +
  geom_point() +
  geom_errorbar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](figures/ldlc-endpoint-1.png)<!-- -->

```r
eval.data %>% 
  summarize_at(.vars=vars(Pct.Wt.Change,`Weight Change`,`Change in LDL-C`,`BMI Change`),
               .funs=funs(mean(.,na.rm=TRUE))) 
```

```
## # A tibble: 1 × 4
##   Pct.Wt.Change `Weight Change` `Change in LDL-C` `BMI Change`
##           <dbl>           <dbl>             <dbl>        <dbl>
## 1         -6.22           -5.76              11.2        -2.24
```

```r
#Above is non-weighted average
```

# Meta-Analysis 


```r
meta.data <-
  eval.data %>%
  filter(!is.na(`Endpoint SD`)) %>%
  mutate(Pooled.SD=sqrt(`Baseline SD`^2+`Endpoint SD`^2)) %>%
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



forest.meta(ldl.c.meta, 
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

forest.meta(bmi.meta, 
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

We evaluated 19 studies for this meta-analysis. Using the meta-analysis method, we found fasting blood LDL-C levels were increased 11.474 mg/dL (95% CI: 1.112 to 21.836) after the ketogenic diet intervention compared to pre-intervention levels, with a significant p-value of 0.03. Across these studies, the I<sup>2</sup> is 0, the p-value for Q is 0.997. This is a highly consistent I^2. 

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
|    0.223|      0|

```r
with(eval.data, cor.test(`Change in LDL-C`,`Baseline Weight`, method="spearman")) %>% tidy %>% kable(caption="Correlation between change in LDL-C and baseline weight")
```



Table: Correlation between change in LDL-C and baseline weight

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|   -0.708|      2271|       0|Spearman's rank correlation rho |two.sided   |

```r
lm(`Change in LDL-C`~`Baseline Weight`+`Percent Male`, data=eval.data) %>% tidy %>% kable(caption="Linear model between change in LDL-C and baseline weight, including gender")
```



Table: Linear model between change in LDL-C and baseline weight, including gender

|term              | estimate| std.error| statistic| p.value|
|:-----------------|--------:|---------:|---------:|-------:|
|(Intercept)       |   69.156|      16.0|     4.322|   0.000|
|`Baseline Weight` |   -0.624|       0.2|    -3.127|   0.006|
|`Percent Male`    |   -4.939|      12.4|    -0.399|   0.695|

```r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight`*`Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term                        | df| sumsq| meansq| statistic| p.value|
|:---------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight`           |  1|  3345|   3345|    15.306|   0.002|
|Sex.Group                   |  2|   257|    128|     0.588|   0.569|
|`Baseline Weight`:Sex.Group |  2|   343|    172|     0.785|   0.475|
|Residuals                   | 14|  3059|    219|        NA|      NA|

```r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  3345|   3345|    15.728|   0.001|
|Sex.Group         |  2|   257|    128|     0.604|   0.559|
|Residuals         | 16|  3403|    213|        NA|      NA|

```r
ldl.weight.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline Weight`, data = eval.data)
ldl.weight.baseline.lm %>% tidy %>% kable
```



|term              | estimate| std.error| statistic| p.value|
|:-----------------|--------:|---------:|---------:|-------:|
|(Intercept)       |   71.122|    14.863|      4.79|   0.000|
|`Baseline Weight` |   -0.667|     0.164|     -4.06|   0.001|

```r
ldl.weightsex.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |   65.910|    16.191|     4.071|   0.001|
|`Baseline Weight`    |   -0.537|     0.216|    -2.488|   0.024|
|Sex.GroupMixed       |   -8.132|     9.305|    -0.874|   0.395|
|Sex.GroupMostly Male |  -13.541|    13.070|    -1.036|   0.316|

```r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  3345|   3345|    15.728|   0.001|
|Sex.Group         |  2|   257|    128|     0.604|   0.559|
|Residuals         | 16|  3403|    213|        NA|      NA|

```r
ldl.bmi.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, data = eval.data)
ldl.bmi.baseline.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    59.87|    14.303|      4.19|   0.002|
|`Baseline BMI` |    -1.61|     0.451|     -3.58|   0.004|

```r
ldl.bmi.baseline.sex.aov <- aov(`Change in LDL-C` ~ `Baseline BMI`*`Sex.Group`, data = eval.data)
ldl.bmi.baseline.sex.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI`           |  1|   968|  967.6|     34.85|   0.000|
|Sex.Group                |  2|   400|  199.9|      7.20|   0.016|
|`Baseline BMI`:Sex.Group |  1|   210|  209.6|      7.55|   0.025|
|Residuals                |  8|   222|   27.8|        NA|      NA|

```r
ldl.bmi.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline BMI` + `Sex.Group`, data = eval.data)
ldl.bmi.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI` |  1|   968|    968|     20.17|   0.002|
|Sex.Group      |  2|   400|    200|      4.17|   0.052|
|Residuals      |  9|   432|     48|        NA|      NA|

```r
ldl.bmi.baseline.mixed.lm<- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mixed"))
ldl.bmi.baseline.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mixed groups")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mixed groups

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   19.752|    13.972|      1.41|   0.200|
|`Baseline BMI` |   -0.475|     0.414|     -1.15|   0.288|

```r
ldl.bmi.baseline.female.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mostly Female"))
ldl.bmi.baseline.female.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mostly Female groups")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mostly Female groups

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    95.03|     21.03|      4.52|   0.139|
|`Baseline BMI` |    -2.61|      0.78|     -3.34|   0.185|

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
|`Baseline BMI`                 |  1| 967.6|  967.6|     13.84|   0.005|
|`Normal weight`                |  1| 115.4|  115.4|      1.65|   0.231|
|`Baseline BMI`:`Normal weight` |  1|  86.8|   86.8|      1.24|   0.294|
|Residuals                      |  9| 629.2|   69.9|        NA|      NA|

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
|(Intercept)    |   -1.632|    23.498|    -0.069|   0.947|
|`Baseline BMI` |    0.145|     0.679|     0.213|   0.837|

Lower baseline BMI was associated with an increased change in LDL-C after consumption of a ketogenic diet (r<sup>2</sup> = 0.538, p-value = 0.004). The association with increased LDL-C was consistent with baseline weight, where a lower baseline weight was associated with an increased change in LDL-C (r<sup>2</sup> = 0.478, p-value = 0.001). 

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
|         0.269|          0|

```r
with(eval.data, cor.test(`Change in LDL-C`,Pct.Wt.Change, method="spearman")) %>% tidy %>% kable(caption="Correlation between change in LDL-C and weight change")
```



Table: Correlation between change in LDL-C and weight change

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|    0.519|       548|   0.024|Spearman's rank correlation rho |two.sided   |

```r
ldl.weightsexinter.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change`*`Sex.Group`, data = eval.data)
ldl.weightsexinter.change.aov %>% tidy %>% kable
```



|term                    | df| sumsq| meansq| statistic| p.value|
|:-----------------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change           |  1|   739|  738.7|     4.517|   0.053|
|Sex.Group               |  2|   451|  225.3|     1.378|   0.287|
|Pct.Wt.Change:Sex.Group |  2|   111|   55.4|     0.339|   0.719|
|Residuals               | 13|  2126|  163.5|        NA|      NA|

```r
ldl.weightsex.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change` + `Sex.Group`, data = eval.data)
ldl.weightsex.change.aov %>% tidy %>% kable
```



|term          | df| sumsq| meansq| statistic| p.value|
|:-------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change |  1|   739|    739|      4.95|   0.042|
|Sex.Group     |  2|   451|    225|      1.51|   0.252|
|Residuals     | 15|  2237|    149|        NA|      NA|

```r
ldl.weight.change.lm <- lm(`Change in LDL-C` ~ `Pct.Wt.Change`, data = eval.data)
ldl.weight.change.lm %>% tidy %>% kable
```



|term          | estimate| std.error| statistic| p.value|
|:-------------|--------:|---------:|---------:|-------:|
|(Intercept)   |    19.02|     5.395|      3.53|   0.003|
|Pct.Wt.Change |     1.59|     0.734|      2.16|   0.045|

```r
ldl.bmisexinter.change.aov <- aov(`Change in LDL-C` ~ `BMI Change`*`Sex.Group`, data = eval.data)
ldl.bmisexinter.change.aov %>% tidy %>% kable
```



|term                   | df| sumsq| meansq| statistic| p.value|
|:----------------------|--:|-----:|------:|---------:|-------:|
|`BMI Change`           |  1|   564|  563.8|     17.05|   0.006|
|Sex.Group              |  2|   654|  327.1|      9.89|   0.013|
|`BMI Change`:Sex.Group |  1|   371|  370.6|     11.21|   0.015|
|Residuals              |  6|   198|   33.1|        NA|      NA|

```r
ldl.bmisex.change.aov <- aov(`Change in LDL-C` ~ `BMI Change` + `Sex.Group`, data = eval.data)
ldl.bmisex.change.aov %>% tidy %>% kable
```



|term         | df| sumsq| meansq| statistic| p.value|
|:------------|--:|-----:|------:|---------:|-------:|
|`BMI Change` |  1|   564|  563.8|      6.93|   0.034|
|Sex.Group    |  2|   654|  327.1|      4.03|   0.069|
|Residuals    |  7|   569|   81.3|        NA|      NA|

```r
ldl.bmi.change.lm <- lm(`Change in LDL-C` ~ `BMI Change`, data = eval.data)
ldl.bmi.change.lm %>% tidy %>% kable
```



|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |    20.14|      6.33|      3.18|   0.011|
|`BMI Change` |     4.79|      2.35|      2.04|   0.072|

```r
ldl.bmi.change.mixed.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mixed"))
ldl.bmi.change.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mixed groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mixed groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |    3.191|      4.13|     0.773|   0.475|
|`BMI Change` |    0.336|      1.27|     0.264|   0.802|

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
|`BMI Change` |     19.0|      8.24|      2.31|   0.260|

```r
ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
ldl.pctbmi.change.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    20.52|     7.018|      2.92|   0.017|
|PCT.BMI.Change |     1.67|     0.903|      1.84|   0.098|

```r
ldl.pctbmisexinter.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change`*`Sex.Group`, data = eval.data)
ldl.pctbmisexinter.change.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|PCT.BMI.Change           |  1|   490|  490.1|     11.08|   0.016|
|Sex.Group                |  2|   727|  363.4|      8.21|   0.019|
|PCT.BMI.Change:Sex.Group |  1|   305|  304.8|      6.89|   0.039|
|Residuals                |  6|   265|   44.2|        NA|      NA|

```r
ldl.pctbmisex.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change` + `Sex.Group`, data = eval.data)
ldl.pctbmisex.change.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|PCT.BMI.Change |  1|   490|  490.1|      6.02|   0.044|
|Sex.Group      |  2|   727|  363.4|      4.46|   0.056|
|Residuals      |  7|   570|   81.5|        NA|      NA|

```r
ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
ldl.pctbmi.change.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    20.52|     7.018|      2.92|   0.017|
|PCT.BMI.Change |     1.67|     0.903|      1.84|   0.098|

```r
ldl.pctbmi.change.mixed.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mixed"))
ldl.pctbmi.change.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mixed groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mixed groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |    3.191|      4.13|     0.773|   0.475|
|`BMI Change` |    0.336|      1.27|     0.264|   0.802|

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
|`BMI Change` |     19.0|      8.24|      2.31|   0.260|
Greater BMI decreases over the study period were associated with a smaller increase in LDL-C after consumption of a ketogenic diet, though this did not reach significance (r<sup>2</sup> = 0.315, p-value = 0.072). The association with the change in LDL-C and decrease in BMI was consistent with weight, with change in weight on LDL-C reaching significance (r<sup>2</sup> = 0.216, p-value = 0.045), where greater decreases in weight were associated with lower increases in LDL-C after consumption of a ketogenic diet. Looking at percent BMI change to account for baseline BMI, greater percent change decreases were associated with a lower increase in LDL-C on a ketogenic diet, though this was not significant (r<sup>2</sup> = 0.274, p-value = 0.098).


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
|    0.733|      0|

```r
with(eval.data, cor.test(`Change in LDL-C`,`Baseline LDL`, method="spearman")) %>% tidy %>% kable(caption="Correlation between baseline and change in LDL-C")
```



Table: Correlation between baseline and change in LDL-C

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|   -0.226|      2172|   0.311|Spearman's rank correlation rho |two.sided   |

```r
ldl.baselinesexinter.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` * `Sex.Group`, data = eval.data)
ldl.baselinesexinter.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline LDL`           |  1|   724|    724|      3.61|   0.076|
|Sex.Group                |  2|  2276|   1138|      5.67|   0.014|
|`Baseline LDL`:Sex.Group |  2|  1039|    519|      2.59|   0.106|
|Residuals                | 16|  3210|    201|        NA|      NA|

```r
ldl.baselinesex.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baselinesex.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline LDL` |  1|   724|    724|      3.07|   0.097|
|Sex.Group      |  2|  2276|   1138|      4.82|   0.021|
|Residuals      | 18|  4249|    236|        NA|      NA|

```r
ldl.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline LDL`, data = eval.data)
ldl.baseline.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   39.884|    19.606|      2.03|   0.055|
|`Baseline LDL` |   -0.255|     0.171|     -1.49|   0.152|

```r
ldl.baselinesex.lm <- lm(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baselinesex.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |   53.173|    17.529|      3.03|   0.007|
|`Baseline LDL`       |   -0.229|     0.151|     -1.51|   0.148|
|Sex.GroupMixed       |  -21.094|     7.703|     -2.74|   0.014|
|Sex.GroupMostly Male |  -27.709|    10.929|     -2.54|   0.021|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline LDL` |  1|   724|    724|      3.07|   0.097|
|Sex.Group      |  2|  2276|   1138|      4.82|   0.021|
|Residuals      | 18|  4249|    236|        NA|      NA|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term      | df| sumsq| meansq| statistic| p.value|
|:---------|--:|-----:|------:|---------:|-------:|
|Sex.Group |  2|  2461|   1231|      4.88|   0.019|
|Residuals | 19|  4787|    252|        NA|      NA|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Percent Male`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Percent Male` |  1|  1284|   1284|      4.31|   0.051|
|Residuals      | 20|  5964|    298|        NA|      NA|

Among individuals, baseline LDL-C was not positively correlated with change in LDL-C after consumption of a ketogenic diet and the relationship was not significant (r<sup>2</sup> = 0.1, p-value = 0.152). 

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
## [1] broom_1.0.5   ggrepel_0.9.4 meta_6.5-0    ggplot2_3.4.4 readxl_1.4.3 
## [6] dplyr_1.1.3   tidyr_1.3.0   knitr_1.44   
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.0    xfun_0.40           bslib_0.5.1        
##  [4] purrr_1.0.2         splines_4.2.2       lattice_0.21-9     
##  [7] colorspace_2.1-0    vctrs_0.6.4         generics_0.1.3     
## [10] htmltools_0.5.6.1   mgcv_1.9-0          yaml_2.3.7         
## [13] utf8_1.2.4          rlang_1.1.1         jquerylib_0.1.4    
## [16] pillar_1.9.0        nloptr_2.0.3        glue_1.6.2         
## [19] withr_2.5.2         lifecycle_1.0.3     stringr_1.5.0      
## [22] munsell_0.5.0       gtable_0.3.4        cellranger_1.1.0   
## [25] evaluate_0.22       labeling_0.4.3      fastmap_1.1.1      
## [28] metafor_4.4-0       fansi_1.0.5         Rcpp_1.0.11        
## [31] backports_1.4.1     scales_1.2.1        cachem_1.0.8       
## [34] jsonlite_1.8.7      farver_2.1.1        lme4_1.1-34        
## [37] digest_0.6.33       stringi_1.7.12      CompQuadForm_1.4.3 
## [40] numDeriv_2016.8-1.1 grid_4.2.2          mathjaxr_1.6-0     
## [43] cli_3.6.1           tools_4.2.2         magrittr_2.0.3     
## [46] sass_0.4.7          tibble_3.2.1        pkgconfig_2.0.3    
## [49] MASS_7.3-60         Matrix_1.5-4.1      xml2_1.3.5         
## [52] minqa_1.2.6         rmarkdown_2.25      rstudioapi_0.15.0  
## [55] R6_2.5.1            boot_1.3-28.1       metadat_1.2-0      
## [58] nlme_3.1-163        compiler_4.2.2
```
