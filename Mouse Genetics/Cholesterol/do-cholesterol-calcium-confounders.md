---
title: "Associations with Diversity Outbred Calcium and Cholesterol Levels"
author: "Dave Bridges"
date: "November 27, 2022"
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

To analyze associations between other phenotypes and cholesterol and calcium to identify potential measured confounders.

# Experimental Details

This analysis uses the complete dataset (F01-F425 and M01-M425). 

# Raw Data


```r
phenotype.filename <- 'Svenson_HFD_DO_phenotype_V12.csv'
```



```r
library(readr) #loads the readr package


phenotype.data <- read_csv(phenotype.filename)
#set phenotypes of zero or  to na
phenotype.data[phenotype.data < 0] <- NA

library(forcats)
cholesterol.data <-
  phenotype.data %>%
  mutate(Diet = fct_recode(as.factor(diet),
                           "NCD"="chow",
                           "HFHS"="hf")) %>%
  select(-diet) %>%
  mutate(chol.avg = rowMeans(select(., starts_with("chol")), 
                             na.rm = TRUE))
```

# Correlation Matrix


```r
library("psych")    
correlation.matrix <- cholesterol.data %>%
  select(-gen,-sex,-Diet) %>%
  select(where(is.numeric)) %>% #removes non numeric covariates
  corr.test(method="spearman", use="pairwise",adjust="BH")

correlation.matrix$p %>%
  as.data.frame() %>%
  select(calcium2) %>%
  filter(calcium2<0.05/162) -> calcium.correlates

correlation.matrix$p %>%
  as.data.frame() %>%
  select(chol.avg) %>%
  filter(chol.avg<0.05/162) -> cholesterol.correlates

overlapping.correlates <- intersect(rownames(calcium.correlates),
                                    rownames(cholesterol.correlates))
```


```r
library(venneuler)

venneuler(c("Calcium"=dim(calcium.correlates)[1]+1,
            "Cholesterol"=dim(cholesterol.correlates)[1]+1,
            "Calcium&Cholesterol"=length(overlapping.correlates)+1)) %>%
  plot
```

![](figures/confounders-venn-1.png)<!-- -->

# Testing Moderation by Correlates



```r
base.model <- lm(chol2 ~ Diet + sex + calcium2, data=cholesterol.data)

confounders.data <-
  data.frame("Measure"=overlapping.correlates) %>%
  mutate(Calcium.mod = NA)

for (measure in overlapping.correlates) {
  formula <- as.formula(paste("chol2 ~ ", measure, " + Diet + sex + calcium2"))
  test.model = lm(formula,data=cholesterol.data)
  confounders.data <- 
    confounders.data %>%
    add_row(Measure=measure,
            Calcium.mod=coefficients(test.model)['calcium2'])
}

confounders.data <-
  confounders.data %>% 
  mutate(Calcium.base = coefficients(base.model)['calcium2']) %>%
  filter(!is.na(Calcium.mod)) 

library(ggplot2)
confounders.data %>%
  mutate(Estimate.change = Calcium.mod-Calcium.base) %>%
  ggplot(aes(y=Estimate.change,x=reorder(Measure,Calcium.mod))) +
  geom_bar(stat="identity") +
  labs(y="Change in Calcium Beta",
       x="Phenotype") +
  theme_classic() +
  theme(text=element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](figures/chol-calcium-effect-changes-1.png)<!-- -->

```r
confounders.data %>%
  arrange(Calcium.mod) %>%
  mutate(Estimate.change = Calcium.mod-Calcium.base) %>%
  kable(caption="Base and Moderated Model Estimates")
```



Table: Base and Moderated Model Estimates

|Measure   | Calcium.mod| Calcium.base| Estimate.change|
|:---------|-----------:|------------:|---------------:|
|hdld2     |       0.168|         12.6|         -12.411|
|chol.avg  |       5.154|         12.6|          -7.425|
|nefa2     |      10.482|         12.6|          -2.098|
|tg2       |      10.783|         12.6|          -1.796|
|fat.mri   |      11.205|         12.6|          -1.374|
|gldh2     |      11.662|         12.6|          -0.918|
|chol1     |      11.664|         12.6|          -0.915|
|hdld1     |      11.757|         12.6|          -0.822|
|bw.24     |      11.885|         12.6|          -0.695|
|bw.25     |      11.903|         12.6|          -0.676|
|bw.23     |      11.929|         12.6|          -0.650|
|perc.fat2 |      12.020|         12.6|          -0.559|
|percfat2  |      12.042|         12.6|          -0.538|
|bw.19     |      12.052|         12.6|          -0.527|
|bw.26     |      12.057|         12.6|          -0.523|
|bw.18     |      12.060|         12.6|          -0.519|
|ftm2      |      12.067|         12.6|          -0.512|
|bw.20     |      12.081|         12.6|          -0.498|
|ttm2      |      12.113|         12.6|          -0.466|
|bw.17     |      12.123|         12.6|          -0.456|
|weight2   |      12.155|         12.6|          -0.424|
|glucose2  |      12.164|         12.6|          -0.416|
|leptin    |      12.166|         12.6|          -0.413|
|bw.15     |      12.181|         12.6|          -0.399|
|t.area2   |      12.189|         12.6|          -0.390|
|bw.22     |      12.199|         12.6|          -0.381|
|bw.21     |      12.200|         12.6|          -0.379|
|bw.16     |      12.217|         12.6|          -0.362|
|bw.5      |      12.291|         12.6|          -0.288|
|percfat1  |      12.366|         12.6|          -0.213|
|perc.fat1 |      12.369|         12.6|          -0.210|
|ftm1      |      12.406|         12.6|          -0.174|
|calcium2  |      12.579|         12.6|           0.000|
|chol2     |      12.579|         12.6|           0.000|

# Potential Mediators

These are defined as phenotypes that associate with both calcium and cholesterol, and *weaken* the sex and diet adjusted association of cholesterol and calcium when included in the model.

## Fat Mass


```r
library(broom)
complete.data <- 
  cholesterol.data %>%
  filter(!(is.na(calcium2))) %>%
  filter(!(is.na(chol2))) %>%
  filter(!(is.na(percfat2)))
         


base.lm <- lm(chol2 ~ Diet + sex + calcium2 + percfat2, data=complete.data)        
fat.lm <- lm(chol2 ~ Diet + sex + calcium2 + percfat2, data=complete.data)
fat.lm %>% tidy %>% kable(caption="Model including fat mass")
```



Table: Model including fat mass

|term        | estimate| std.error| statistic| p.value|
|:-----------|--------:|---------:|---------:|-------:|
|(Intercept) |    -43.0|     8.074|     -5.32|       0|
|DietHFHS    |     22.3|     2.111|     10.59|       0|
|sexM        |     19.8|     1.805|     10.97|       0|
|calcium2    |     12.0|     0.866|     13.90|       0|
|percfat2    |     59.0|    11.328|      5.21|       0|

```r
library(mediation)
results <- mediate(fat.lm, base.lm, treat='calcium2', mediator='percfat2',
                   boot=TRUE, sims=1000)
results %>% summary
```

```
## 
## Causal Mediation Analysis 
## 
## Nonparametric Bootstrap Confidence Intervals with the Percentile Method
## 
##                Estimate 95% CI Lower 95% CI Upper p-value    
## ACME            710.166      454.995       964.55  <2e-16 ***
## ADE              12.042       10.404        13.61  <2e-16 ***
## Total Effect    722.208      467.042       978.03  <2e-16 ***
## Prop. Mediated    0.983        0.974         0.99  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Sample Size Used: 758 
## 
## 
## Simulations: 1000
```


# Potential Confounders

These are defined as phenotypes that associate with both calcium and cholesterol, and *strengthen* the sex and diet adjusted association of cholesterol and calcium when included in the model.

## Phosphorous


```r
phos.lm <- lm(chol2 ~ Diet + sex + calcium2 + phosphorus2, data=cholesterol.data)
library(broom)
phos.lm %>% tidy %>% kable(caption="Model including phosphorus")
```



Table: Model including phosphorus

|term        | estimate| std.error| statistic| p.value|
|:-----------|--------:|---------:|---------:|-------:|
|(Intercept) |   -31.89|     8.003|     -3.98|   0.000|
|DietHFHS    |    28.05|     1.801|     15.57|   0.000|
|sexM        |    18.86|     1.839|     10.26|   0.000|
|calcium2    |    13.64|     0.972|     14.04|   0.000|
|phosphorus2 |    -2.07|     0.854|     -2.42|   0.016|

```r
#anova(base.model,phos.lm) %>% tidy %>% kable(caption="Chi squared test of models with or without calcium")
```

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
##  [1] mediation_4.5.0 sandwich_3.0-2  mvtnorm_1.1-3   Matrix_1.5-3   
##  [5] MASS_7.3-58.1   broom_1.0.2     ggplot2_3.4.0   venneuler_1.1-3
##  [9] rJava_1.0-6     psych_2.2.9     forcats_0.5.2   readr_2.1.3    
## [13] dplyr_1.0.10    tidyr_1.2.1     knitr_1.41     
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-161        bit64_4.0.5         RColorBrewer_1.1-3 
##  [4] tools_4.2.2         backports_1.4.1     bslib_0.4.2        
##  [7] utf8_1.2.2          R6_2.5.1            rpart_4.1.19       
## [10] Hmisc_4.7-2         DBI_1.1.3           colorspace_2.0-3   
## [13] nnet_7.3-18         withr_2.5.0         tidyselect_1.2.0   
## [16] gridExtra_2.3       mnormt_2.1.1        bit_4.0.5          
## [19] compiler_4.2.2      cli_3.6.0           htmlTable_2.4.1    
## [22] labeling_0.4.2      sass_0.4.4          scales_1.2.1       
## [25] checkmate_2.1.0     stringr_1.5.0       digest_0.6.31      
## [28] foreign_0.8-84      minqa_1.2.5         rmarkdown_2.19     
## [31] base64enc_0.1-3     jpeg_0.1-10         pkgconfig_2.0.3    
## [34] htmltools_0.5.4     lme4_1.1-31         fastmap_1.1.0      
## [37] highr_0.10          htmlwidgets_1.6.1   rlang_1.0.6        
## [40] rstudioapi_0.14     jquerylib_0.1.4     farver_2.1.1       
## [43] generics_0.1.3      zoo_1.8-11          jsonlite_1.8.4     
## [46] vroom_1.6.0         magrittr_2.0.3      Formula_1.2-4      
## [49] interp_1.1-3        Rcpp_1.0.9          munsell_0.5.0      
## [52] fansi_1.0.3         lifecycle_1.0.3     stringi_1.7.12     
## [55] yaml_2.3.6          grid_4.2.2          parallel_4.2.2     
## [58] crayon_1.5.2        deldir_1.0-6        lattice_0.20-45    
## [61] splines_4.2.2       hms_1.1.2           pillar_1.8.1       
## [64] boot_1.3-28.1       lpSolve_5.6.17      glue_1.6.2         
## [67] evaluate_0.19       latticeExtra_0.6-30 data.table_1.14.6  
## [70] nloptr_2.0.3        png_0.1-8           vctrs_0.5.1        
## [73] tzdb_0.3.0          gtable_0.3.1        purrr_1.0.1        
## [76] assertthat_0.2.1    cachem_1.0.6        xfun_0.36          
## [79] survival_3.5-0      tibble_3.1.8        cluster_2.1.4      
## [82] ellipsis_0.3.2
```

