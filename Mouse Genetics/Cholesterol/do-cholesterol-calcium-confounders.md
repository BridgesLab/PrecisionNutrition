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
phenotype.filename <- 'Svenson-183_Svenson_DO-phenotypes.csv'
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
  select(-coat_color,-gen,-sex,-Diet,-sample) %>%
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

|Measure     | Calcium.mod| Calcium.base| Estimate.change|
|:-----------|-----------:|------------:|---------------:|
|hdld2       |       0.424|         12.7|         -12.236|
|chol.avg    |       5.276|         12.7|          -7.383|
|nefa2       |      10.411|         12.7|          -2.249|
|tg2         |      10.801|         12.7|          -1.858|
|fat_mri     |      11.205|         12.7|          -1.454|
|gldh2       |      11.713|         12.7|          -0.947|
|chol1       |      11.840|         12.7|          -0.819|
|bw_24       |      11.911|         12.7|          -0.749|
|hdld1       |      11.921|         12.7|          -0.738|
|bw_25       |      11.927|         12.7|          -0.733|
|bw_23       |      11.952|         12.7|          -0.708|
|bw_26       |      12.082|         12.7|          -0.578|
|bw_19       |      12.158|         12.7|          -0.502|
|bw_18       |      12.166|         12.7|          -0.494|
|perc_fat2   |      12.178|         12.7|          -0.482|
|bw_20       |      12.185|         12.7|          -0.475|
|ftm2        |      12.211|         12.7|          -0.448|
|bw_22       |      12.221|         12.7|          -0.439|
|bw_17       |      12.227|         12.7|          -0.432|
|glucose2    |      12.232|         12.7|          -0.427|
|leptin      |      12.262|         12.7|          -0.398|
|ttm2        |      12.264|         12.7|          -0.396|
|bw_15       |      12.291|         12.7|          -0.368|
|weight2     |      12.306|         12.7|          -0.354|
|bw_21       |      12.339|         12.7|          -0.321|
|t_area2     |      12.341|         12.7|          -0.319|
|bw_5        |      12.350|         12.7|          -0.310|
|perc_fat1   |      12.451|         12.7|          -0.209|
|ftm1        |      12.498|         12.7|          -0.162|
|chol2       |      12.660|         12.7|           0.000|
|calcium2    |      12.660|         12.7|           0.000|
|phosphorus2 |      13.705|         12.7|           1.045|

# Potential Mediators

These are defined as phenotypes that associate with both calcium and cholesterol, and *weaken* the sex and diet adjusted association of cholesterol and calcium when included in the model.

## Fat Mass


```r
library(broom)
complete.data <- 
  cholesterol.data %>%
  filter(!(is.na(calcium2))) %>%
  filter(!(is.na(chol2))) %>%
  filter(!(is.na(fat_mri)))
         


base.lm <- lm(chol2 ~ Diet + sex + calcium2 + fat_mri, data=complete.data)        
fat.lm <- lm(chol2 ~ Diet + sex + calcium2 + fat_mri, data=complete.data)
fat.lm %>% tidy %>% kable(caption="Model including fat mass")
```



Table: Model including fat mass

|term        | estimate| std.error| statistic| p.value|
|:-----------|--------:|---------:|---------:|-------:|
|(Intercept) |  -22.978|    21.272|    -1.080|   0.282|
|DietHFHS    |   35.744|     5.313|     6.728|   0.000|
|sexM        |   10.953|     4.662|     2.349|   0.020|
|calcium2    |   11.205|     2.382|     4.705|   0.000|
|fat_mri     |    0.218|     0.536|     0.407|   0.685|

```r
library(mediation)
results <- mediate(fat.lm, base.lm, treat='calcium2', mediator='fat_mri',
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
## ACME              2.447      -14.722        13.93   0.740    
## ADE              11.205        6.836        16.14  <2e-16 ***
## Total Effect     13.652       -1.844        25.66   0.076 .  
## Prop. Mediated    0.179       -5.480         4.45   0.664    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Sample Size Used: 148 
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
|(Intercept) |   -32.45|     7.975|     -4.07|   0.000|
|DietHFHS    |    27.95|     1.797|     15.55|   0.000|
|sexM        |    19.03|     1.836|     10.36|   0.000|
|calcium2    |    13.71|     0.962|     14.25|   0.000|
|phosphorus2 |    -2.09|     0.849|     -2.46|   0.014|

```r
#anova(base.model,phos.lm) %>% tidy %>% kable(caption="Chi squared test of models with or without calcium")
```

# Session Information


```r
sessionInfo()
```

```
## R version 4.2.0 (2022-04-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur/Monterey 10.16
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
##  [5] MASS_7.3-58.1   broom_1.0.1     ggplot2_3.4.0   venneuler_1.1-3
##  [9] rJava_1.0-6     psych_2.2.9     forcats_0.5.2   readr_2.1.3    
## [13] dplyr_1.0.10    tidyr_1.2.1     knitr_1.41     
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-160        bit64_4.0.5         RColorBrewer_1.1-3 
##  [4] tools_4.2.0         backports_1.4.1     bslib_0.4.1        
##  [7] utf8_1.2.2          R6_2.5.1            rpart_4.1.19       
## [10] Hmisc_4.7-2         DBI_1.1.3           colorspace_2.0-3   
## [13] nnet_7.3-18         withr_2.5.0         tidyselect_1.2.0   
## [16] gridExtra_2.3       mnormt_2.1.1        bit_4.0.5          
## [19] compiler_4.2.0      cli_3.4.1           htmlTable_2.4.1    
## [22] labeling_0.4.2      sass_0.4.3          scales_1.2.1       
## [25] checkmate_2.1.0     stringr_1.4.1       digest_0.6.30      
## [28] foreign_0.8-83      minqa_1.2.5         rmarkdown_2.18     
## [31] base64enc_0.1-3     jpeg_0.1-9          pkgconfig_2.0.3    
## [34] htmltools_0.5.3     lme4_1.1-31         fastmap_1.1.0      
## [37] highr_0.9           htmlwidgets_1.5.4   rlang_1.0.6        
## [40] rstudioapi_0.14     jquerylib_0.1.4     farver_2.1.1       
## [43] generics_0.1.3      zoo_1.8-11          jsonlite_1.8.3     
## [46] vroom_1.6.0         magrittr_2.0.3      Formula_1.2-4      
## [49] interp_1.1-3        Rcpp_1.0.9          munsell_0.5.0      
## [52] fansi_1.0.3         lifecycle_1.0.3     stringi_1.7.8      
## [55] yaml_2.3.6          grid_4.2.0          parallel_4.2.0     
## [58] crayon_1.5.2        deldir_1.0-6        lattice_0.20-45    
## [61] splines_4.2.0       hms_1.1.2           pillar_1.8.1       
## [64] boot_1.3-28         lpSolve_5.6.17      glue_1.6.2         
## [67] evaluate_0.18       latticeExtra_0.6-30 data.table_1.14.6  
## [70] nloptr_2.0.3        png_0.1-7           vctrs_0.5.1        
## [73] tzdb_0.3.0          gtable_0.3.1        purrr_0.3.5        
## [76] assertthat_0.2.1    cachem_1.0.6        xfun_0.35          
## [79] survival_3.4-0      tibble_3.1.8        cluster_2.1.4      
## [82] ellipsis_0.3.2
```

