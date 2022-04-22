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

Reviewed data from the Choi *et al* meta-analysis (http://dx.doi.org/10.3390/nu12072005), pulling in data on baseline weight, weight changes, LDL, LDL changes and standard deviations. A systematic literature search of PubMed was then performed to identify other randomized controlled trials (RCTs) and single-arm interventions of patients that evaluated the effects of a ketogenic diet on weight and lipid profile as primary endpoints. All studies using a KD diet that met our inclusion criteria where intake of carbohydrate was less than 25 grams per day were included. This search was most recently updated on Fri Apr 22 09:16:03 2022.

We used a value 130mg/dL of LDL-C at baseline to stratify individuals as being hypercholesterolemic or not.


Correlations in this meta-analysis before and after the administration of the ketogenic diet were analyzed with linear models and results given using pearsons correlation coefficient, statistical significance was defined as below 0.05. 


For all outcomes, we tested sex as a modifier and as a covariate. For outcomes where sex was found to be a significant modifier, these results are reported. 


```r
library(readxl) #loads the readr package
filename <- 'LDL Study Summary.xlsx' #make this a separate line, you can use any variable you want

#this loads whatever the file is into a dataframe called exp.data if it exists
exp.data <- read_excel(filename)
eval.data <- exp.data %>% 
  filter(Use=='x') %>%
  mutate(Pct.Wt.Change = `Weight Change`/`Baseline Weight`*100)%>%
  mutate(PCT.BMI.Change = `BMI Change` / `Baseline BMI`*100)%>%
  mutate(Sex.Group = cut(`Percent Male`, breaks = c(0,.1,.9,1), include.lowest = TRUE, labels = c("Mostly Female", "Mixed", "Mostly Male")))
```

These data can be found in **C:/Users/Cody/OneDrive/Documents/GitHub/PrecisionNutrition/Meta Analysis** in a file named **LDL Study Summary.xlsx**.  This script was most recently updated on **Fri Apr 22 09:16:04 2022**.

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
## # A tibble: 1 x 4
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
  mutate(SMD.Wt = `Weight Change`) 

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

#forest.meta(ldl.c.meta, 
            #sortvar = SMD,
            #predict = F, 
            #print.tau2 = TRUE,
           # print.I2 = TRUE,
            #leftcols=c('Study'),
            #rightlabs=c('Change','95% CI','Weight'),
           # layout = "RevMan5",
          #  fontsize=8,
          #  leftcols = c("Study",'effect',"seTE",'ci','w.fixed'),
          #  leftlabs = c("Study","Change","SE","95% CI","Weight"),
          #  range = F,
          #  digits=1,
          #  digits.se=1)
```

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
|(Intercept)       |   69.286|    16.562|     4.183|   0.001|
|`Baseline Weight` |   -0.627|     0.208|    -3.013|   0.008|
|`Percent Male`    |   -4.849|    12.797|    -0.379|   0.710|

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
|(Intercept)          |   65.991|    16.716|     3.948|   0.001|
|`Baseline Weight`    |   -0.538|     0.223|    -2.415|   0.029|
|Sex.GroupMixed       |   -8.305|     9.667|    -0.859|   0.404|
|Sex.GroupMostly Male |  -13.512|    13.489|    -1.002|   0.332|

```r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  3319|   3319|    14.655|   0.002|
|Sex.Group         |  2|   259|    130|     0.572|   0.576|
|Residuals         | 15|  3397|    226|        NA|      NA|

```r
ldl.bmi.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, data = eval.data)
ldl.bmi.baseline.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    59.87|    14.303|      4.19|   0.002|
|`Baseline BMI` |    -1.61|     0.451|     -3.58|   0.004|

```r
ldl.bmi.baseline.sex.aov <- aov(`Change in LDL-C` ~ `Baseline BMI` + `Sex.Group`, data = eval.data)
ldl.bmi.baseline.sex.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI` |  1|   965|  965.0|     18.51|   0.003|
|Sex.Group      |  2|   411|  205.3|      3.94|   0.065|
|Residuals      |  8|   417|   52.1|        NA|      NA|

```r
ldl.weight.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weight.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  3319|   3319|    14.655|   0.002|
|Sex.Group         |  2|   259|    130|     0.572|   0.576|
|Residuals         | 15|  3397|    226|        NA|      NA|

```r
ldl.bmi.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline BMI` + `Sex.Group`, data = eval.data)
ldl.bmi.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI` |  1|   965|  965.0|     18.51|   0.003|
|Sex.Group      |  2|   411|  205.3|      3.94|   0.065|
|Residuals      |  8|   417|   52.1|        NA|      NA|

```r
ldl.bmi2.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline BMI`*`Sex.Group`, data = eval.data)
ldl.bmi2.baseline.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI`           |  1|   965|  965.0|     32.26|   0.001|
|Sex.Group                |  2|   411|  205.3|      6.86|   0.022|
|`Baseline BMI`:Sex.Group |  1|   208|  207.7|      6.94|   0.034|
|Residuals                |  7|   209|   29.9|        NA|      NA|

```r
ldl.weight.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight`, data = eval.data)
ldl.weight.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  3345|   3345|      16.5|   0.001|
|Residuals         | 18|  3660|    203|        NA|      NA|

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
             x=`PCT.BMI.Change`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="% BMI Change")
```

![](figures/ldl-change-vs-weight-change-3.png)<!-- -->

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

![](figures/ldl-change-vs-weight-change-4.png)<!-- -->

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
ldl.weight.change.lm <- lm(`Change in LDL-C` ~ `Pct.Wt.Change`, data = eval.data)
ldl.weight.change.lm %>% tidy %>% kable
```



|term          | estimate| std.error| statistic| p.value|
|:-------------|--------:|---------:|---------:|-------:|
|(Intercept)   |    19.02|     5.395|      3.53|   0.003|
|Pct.Wt.Change |     1.59|     0.734|      2.16|   0.045|

```r
ldl.weightsex.change.lm <- lm(`Change in LDL-C` ~ `Pct.Wt.Change` + `Sex.Group`, data = eval.data)
ldl.weightsex.change.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |    24.22|     6.314|      3.84|   0.002|
|Pct.Wt.Change        |     1.31|     0.878|      1.49|   0.158|
|Sex.GroupMixed       |    -8.02|     7.977|     -1.00|   0.332|
|Sex.GroupMostly Male |   -17.37|    10.604|     -1.64|   0.124|

```r
ldl.weightsex.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change` + `Sex.Group`, data = eval.data)
ldl.weightsex.change.aov %>% tidy %>% kable
```



|term          | df| sumsq| meansq| statistic| p.value|
|:-------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change |  1|   735|    735|      4.60|   0.050|
|Sex.Group     |  2|   449|    225|      1.40|   0.278|
|Residuals     | 14|  2237|    160|        NA|      NA|

```r
ldl.bmi.change.lm <- lm(`Change in LDL-C` ~ `BMI Change`, data = eval.data)
ldl.bmi.change.lm %>% tidy %>% kable
```



|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |    20.14|      6.33|      3.18|   0.011|
|`BMI Change` |     4.79|      2.35|      2.04|   0.072|

```r
ldl.bmisex.change.lm <- lm(`Change in LDL-C` ~ `BMI Change` + `Sex.Group`, data = eval.data)
ldl.bmisex.change.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |    27.67|      5.76|     4.801|   0.002|
|`BMI Change`         |     1.66|      2.25|     0.739|   0.484|
|Sex.GroupMixed       |   -20.64|      7.42|    -2.781|   0.027|
|Sex.GroupMostly Male |   -15.33|     10.41|    -1.472|   0.185|

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
ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
ldl.pctbmi.change.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    20.52|     7.018|      2.92|   0.017|
|PCT.BMI.Change |     1.67|     0.903|      1.84|   0.098|

```r
ldl.pctbmisex.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change` + `Sex.Group`, data = eval.data)
ldl.pctbmisex.change.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |   28.078|     6.051|     4.640|   0.002|
|PCT.BMI.Change       |    0.578|     0.794|     0.728|   0.490|
|Sex.GroupMixed       |  -21.108|     7.123|    -2.963|   0.021|
|Sex.GroupMostly Male |  -15.077|    10.423|    -1.447|   0.191|

```r
ldl.weight.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change` + `Sex.Group`, data = eval.data)
ldl.weight.change.aov %>% tidy %>% kable
```



|term          | df| sumsq| meansq| statistic| p.value|
|:-------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change |  1|   735|    735|      4.60|   0.050|
|Sex.Group     |  2|   449|    225|      1.40|   0.278|
|Residuals     | 14|  2237|    160|        NA|      NA|

```r
ldl.bmi.change.aov <- aov(`Change in LDL-C` ~ `BMI Change`, data = eval.data)
ldl.bmi.change.aov %>% tidy %>% kable
```



|term         | df| sumsq| meansq| statistic| p.value|
|:------------|--:|-----:|------:|---------:|-------:|
|`BMI Change` |  1|   564|    564|      4.15|   0.072|
|Residuals    |  9|  1223|    136|        NA|      NA|

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
ldl.pctbmi.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change` + `Sex.Group`, data = eval.data)
ldl.pctbmi.change.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|PCT.BMI.Change |  1|   490|  490.1|      6.02|   0.044|
|Sex.Group      |  2|   727|  363.4|      4.46|   0.056|
|Residuals      |  7|   570|   81.5|        NA|      NA|

```r
ldl.weight.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change`, data = eval.data)
ldl.weight.change.aov %>% tidy %>% kable
```



|term          | df| sumsq| meansq| statistic| p.value|
|:-------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change |  1|   739|    739|      4.67|   0.045|
|Residuals     | 17|  2687|    158|        NA|      NA|
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
|(Intercept)          |   53.821|    18.117|      2.97|   0.009|
|`Baseline LDL`       |   -0.235|     0.157|     -1.50|   0.152|
|Sex.GroupMixed       |  -21.424|     7.980|     -2.68|   0.016|
|Sex.GroupMostly Male |  -27.756|    11.217|     -2.47|   0.024|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline LDL` |  1|   706|    706|      2.84|   0.110|
|Sex.Group      |  2|  2298|   1149|      4.62|   0.025|
|Residuals      | 17|  4226|    249|        NA|      NA|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term      | df| sumsq| meansq| statistic| p.value|
|:---------|--:|-----:|------:|---------:|-------:|
|Sex.Group |  2|  2446|   1223|       4.6|   0.024|
|Residuals | 18|  4784|    266|        NA|      NA|

```r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Percent Male`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Percent Male` |  1|  1293|   1293|      4.14|   0.056|
|Residuals      | 19|  5937|    312|        NA|      NA|

Among individuals, baseline LDL-C was not positively correlated with change in LDL-C after consumption of a ketogenic diet and the relationship was not significant (r<sup>2</sup> = 0.1, p-value = 0.152). 

# Session Information


```r
sessionInfo()
```

```
## R version 4.1.3 (2022-03-10)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 22000)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] broom_0.8.0   ggrepel_0.9.1 meta_5.2-0    ggplot2_3.3.5 readxl_1.4.0 
## [6] dplyr_1.0.8   tidyr_1.2.0   knitr_1.38   
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.2   xfun_0.30          bslib_0.3.1        purrr_0.3.4       
##  [5] splines_4.1.3      lattice_0.20-45    colorspace_2.0-3   vctrs_0.4.1       
##  [9] generics_0.1.2     htmltools_0.5.2    mgcv_1.8-39        yaml_2.3.5        
## [13] utf8_1.2.2         rlang_1.0.2        jquerylib_0.1.4    pillar_1.7.0      
## [17] nloptr_2.0.0       glue_1.6.2         withr_2.5.0        lifecycle_1.0.1   
## [21] stringr_1.4.0      munsell_0.5.0      gtable_0.3.0       cellranger_1.1.0  
## [25] evaluate_0.15      labeling_0.4.2     fastmap_1.1.0      metafor_3.4-0     
## [29] fansi_1.0.3        highr_0.9          Rcpp_1.0.8.3       backports_1.4.1   
## [33] scales_1.2.0       jsonlite_1.8.0     farver_2.1.0       lme4_1.1-29       
## [37] digest_0.6.29      stringi_1.7.6      CompQuadForm_1.4.3 mathjaxr_1.6-0    
## [41] grid_4.1.3         cli_3.2.0          tools_4.1.3        magrittr_2.0.3    
## [45] sass_0.4.1         tibble_3.1.6       crayon_1.5.1       pkgconfig_2.0.3   
## [49] ellipsis_0.3.2     MASS_7.3-55        Matrix_1.4-0       xml2_1.3.3        
## [53] minqa_1.2.4        rmarkdown_2.13     rstudioapi_0.13    metadat_1.2-0     
## [57] R6_2.5.1           boot_1.3-28        nlme_3.1-155       compiler_4.1.3
```
