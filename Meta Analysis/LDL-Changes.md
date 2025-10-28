---
title: "Systematic Review of LDL Changes on a Ketogenic Diet"
author: "Cody Cousineau, Kai-Lin Jen and Dave Bridges"
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

<<<<<<< Updated upstream
Reviewed data from the Choi *et al* meta-analysis (http://dx.doi.org/10.3390/nu12072005), pulling in data on baseline weight, weight changes, LDL, LDL changes and standard deviations. A systematic literature search of PubMed was then performed to identify other randomized controlled trials (RCTs) and single-arm interventions of patients that evaluated the effects of a ketogenic diet on weight and lipid profile as primary endpoints. All studies using a KD diet that met our inclusion criteria where intake of carbohydrate was less than 25 grams per day were included. This search was most recently updated on Fri Aug  2 15:48:31 2024.
=======
Reviewed data from the Choi *et al* meta-analysis (http://dx.doi.org/10.3390/nu12072005), pulling in data on baseline weight, weight changes, LDL, LDL changes and standard deviations. A systematic literature search of PubMed was then performed to identify other randomized controlled trials (RCTs) and single-arm interventions of patients that evaluated the effects of a ketogenic diet on weight and lipid profile as primary endpoints. All studies using a KD diet that met our inclusion criteria where intake of carbohydrate was less than 25 grams per day were included. This search was most recently updated on Tue Oct 28 09:43:20 2025.
>>>>>>> Stashed changes

We used a value 130mg/dL of LDL-C at baseline to stratify individuals as being hypercholesterolemic or not.


Correlations in this meta-analysis before and after the administration of the ketogenic diet were analyzed with linear models and results given using pearsons correlation coefficient, statistical significance was defined as below 0.05. 


For all outcomes, we tested sex as a modifier and as a covariate. For outcomes where sex was found to be a significant modifier, these results are reported. 


``` r
<<<<<<< Updated upstream
filename <- 'VLCF Meta-Analysis.csv'
#filename <- 'LDL Study Summary.xlsx' #make this a separate line, you can use any variable you want
google.sheet.link <- 'https://docs.google.com/spreadsheets/d/e/2PACX-1vRkVL0hHHPMJ_fI8EVc_RiFJvcnTHmBySUd3bvikUJjancG8DCsaG5k0eIFUBfov4drTr8MMNV8GLnv/pub?gid=0&single=true&output=csv'
download.file(google.sheet.link,destfile=filename) #can commnt this out if you dont want to udpat the file
=======
library(readr) #loads the readr package
filename <- 'Ketogenic Diet Meta-Analysis Data - Sheet1.csv' #make this a separate line, you can use any variable you want
>>>>>>> Stashed changes

#this loads whatever the file is into a dataframe called exp.data if it exists
#make a new variable called break.groups, make it so c is the break.groups so that is the code

<<<<<<< Updated upstream
exp.data <- read_csv(filename)
eval.data <- 
  exp.data %>% 
=======
exp.data <- read_csv(filename,col_types=cols(
  Study = col_character(),
  `Baseline Weight` = col_double(),
  `Baseline Weight SD` = col_double(),
  `Baseline Weight CI` = col_double(),
  `Weight Change` = col_double(),
  `Weight Change SD` = col_double(),
  `Weight Change CI` = col_double(),
  `Baseline BMI` = col_double(),
  `Baseline BMI SD` = col_double(),
  `Baseline BMI CI` = col_double(),
  `BMI Change` = col_double(),
  `BMI Change SD` = col_double(),
  `BMI Change CI` = col_double(),
  `Baseline LDL` = col_double(),
  `Baseline LDL SD` = col_double(),
  `Baseline LDL CI` = col_double(),
  `Endpoint LDL` = col_double(),
  `LDL Endpoint SD` = col_double(),
  `LDL Endpoint CI` = col_double(),
  `Change in LDL-C` = col_double(),
  `Change in LDL SD` = col_double(),
  `Change in LDL CI` = col_double(),
  `Baseline TG` = col_double(),
  `Baseline TG SD` = col_double(),
  `Baseline TG CI` = col_double(),
  `Endpoint TG` = col_double(),
  `Endpoint TG SD` = col_double(),
  `Endpoint TG CI` = col_double(),
  `Change in TG` = col_double(),
  `Change in TG SD` = col_double(),
  `Change in TG CI` = col_double(),
  `Baseline HOMA-IR` = col_double(),
  `Baseline HOMA-IR SD` = col_double(),
  `Baseline HOMA-IR CI` = col_double(),
  `Endpoint HOMA-IR` = col_double(),
  `Endpoint HOMA-IR SD` = col_double(),
  `Endpoint HOMA-IR CI` = col_double(),
  `Baseline SBP` = col_double(),
  `Baseline SBP SD` = col_double(),
  `Baseline SBP CI` = col_double(),
  `Endpoint SBP` = col_double(),
  `Endpoint SBP SD` = col_double(),
  `Endpoint SBP CI` = col_double(),
  `change in SBP` = col_double(),
  `change in SBP SD` = col_double(),
  `change in SBP CI` = col_double(),
  `Baseline DBP` = col_double(),
  `Baseline DBP SD` = col_double(),
  `Baseline DBP CI` = col_double(),
  `change in DBP` = col_double(),
  `change in DBP SD` = col_double(),
  `change in DBP CI` = col_double(),
  `Endpoint DBP` = col_double(),
  `Endpoint DBP SD` = col_double(),
  `Endpoint DBP CI` = col_double(),
  n = col_double(),
  Time = col_double(),
  `Inter CI` = col_double(),
  `Inter SD` = col_double(),
  Notes = col_character(),
  Use = col_character(),
  `Percent Male` = col_double(),
  `Normal weight` = col_character(),
  Link = col_character()
))
eval.data <- exp.data %>% 
>>>>>>> Stashed changes
  filter(Use=='x') %>%
  mutate(across(`Baseline Weight`:n,.fns=as.numeric)) %>% #force all fields to be numeric
  mutate(Pct.Wt.Change = `Weight Change`/`Baseline Weight`*100)%>%
  mutate(PCT.BMI.Change = `BMI Change` / `Baseline BMI`*100)%>%
  mutate(Sex.Group = cut(`Percent Male`, breaks = c(0,.1,.9,1), include.lowest = TRUE, labels = c("Mostly Female", "Mixed", "Mostly Male")))
```

<<<<<<< Updated upstream
These data can be found in **/Users/davebrid/Documents/GitHub/PrecisionNutrition/Meta Analysis** in a file named **VLCF Meta-Analysis.csv**.  This script was most recently updated on **Fri Aug  2 15:48:34 2024**.
=======
These data can be found in **/Users/davebrid/Documents/GitHub/PrecisionNutrition/Meta Analysis** in a file named **Ketogenic Diet Meta-Analysis Data - Sheet1.csv**.  This script was most recently updated on **Tue Oct 28 09:43:21 2025**.
>>>>>>> Stashed changes

This analysis includes 22 studies with 463 total participants.


<<<<<<< Updated upstream
=======
``` r
library(ggplot2)

eval.data <-
  eval.data %>%
  mutate(LDL.endpoint = `Baseline LDL`+`Change in LDL-C`)

eval.data %>%
  ggplot(aes(y=LDL.endpoint,
             ymin=LDL.endpoint-`LDL Endpoint SD`,
             ymax=LDL.endpoint+`LDL Endpoint SD`,
             x=reorder(Study,-`LDL.endpoint` ))) +
  geom_point() +
  geom_errorbar()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](figures/ldlc-endpoint-1.png)<!-- -->

``` r
eval.data %>% 
  summarize_at(.vars=vars(Pct.Wt.Change,`Weight Change`,`Change in LDL-C`,`BMI Change`),
               .funs=funs(mean(.,na.rm=TRUE))) 
```

```
## # A tibble: 1 × 4
##   Pct.Wt.Change `Weight Change` `Change in LDL-C` `BMI Change`
##           <dbl>           <dbl>             <dbl>        <dbl>
## 1         -6.33           -5.95              8.00        -2.51
```

``` r
#Above is non-weighted average
```
>>>>>>> Stashed changes

# Meta-Analysis 


``` r
<<<<<<< Updated upstream
meta.data <-
  eval.data %>%
  filter(!is.na(`LDL Endpoint SD`)) %>%
  mutate(Pooled.LDL.SD=sqrt(`Baseline LDL SD`^2+`LDL Endpoint SD`^2)) %>%
  mutate(SMD=`Change in LDL-C`) %>%
=======
ldl.meta.data <-
  eval.data %>%
  mutate(`Baseline LDL SD` = case_when(is.na(`Baseline LDL SD`)~`Baseline LDL CI`/1.96,
                                      TRUE~`Baseline LDL SD`)) %>%
    mutate(`LDL Endpoint SD` = case_when(is.na(`LDL Endpoint SD`)~`LDL Endpoint CI`/1.96,
                                      TRUE~`LDL Endpoint SD`)) %>%
  mutate(Pooled.SD.LDL=sqrt(`Baseline LDL SD`^2+`LDL Endpoint SD`^2)) %>%
  filter(!is.na(Pooled.SD.LDL)) %>%
  #filter(is.na(`Change in LDL-C`)) %>%
  mutate(MD.LDL=`Change in LDL-C`) %>%
>>>>>>> Stashed changes
  mutate(SMD.Wt = `Weight Change`) %>%
  mutate(Pooled.Wt.SD=sqrt(`Baseline Weight SD`)) 

library(meta)
<<<<<<< Updated upstream
ldl.c.meta <- metagen(TE = `Change in LDL-C`,
                 seTE = Pooled.LDL.SD,
=======
ldl.c.meta <- metagen(TE = MD.LDL,
                 seTE = Pooled.SD.LDL,
>>>>>>> Stashed changes
                 n.e=n,
                 n.c=n,
                 studlab = Study,
                 data = ldl.meta.data,
                 sm = "MD",
                 common=FALSE,
                 random=TRUE,
                 #method.tau = "REML",
                 method.random.ci="HK",
                 title = "LDL-C Changes in Ketogenic Diet Studies")

<<<<<<< Updated upstream
plot(ldl.c.meta)
=======


forest(ldl.c.meta, 
            sortvar = MD.LDL,
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
>>>>>>> Stashed changes
```

![](figures/ldl-meta-analysis-1.png)<!-- -->

<<<<<<< Updated upstream
We evaluated 11 studies for this meta-analysis. Using the meta-analysis method, we found fasting blood LDL-C levels were increased 4.803 mg/dL (95% CI: -20.712 to 30.317) after the ketogenic diet intervention compared to pre-intervention levels, with a significant p-value of 0.712. Across these studies, the I<sup>2</sup> is 0, the p-value for Q is 1. This is a highly consistent I^2. 

## Baysian Meta-Analysis

I wanted to try a Bayesian approach to this, so followed the examples here https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/bayesian-ma.html

The prior probabilities were set using these assumptions

$$\hat\theta_k \sim N(\theta_k, \sigma^2_k)$$
$$\theta_k \sim N(\mu, \tau^2) $$
Presuming a normal distribution with a SD of 44 (the pooled SD but a mean of zero)

$$\mu \sim N(0,44) $$
$$\tau \sim HC(0,0.5)$$


``` r
#library(extraDistr)
#phcauchy(0.3, sigma = 0.3)
library(brms)
mean.sd <- mean(meta.data$Pooled.LDL.SD,na.rm=T)
priors <- c(prior(normal(20,5), class = Intercept), #20 mg incresae seems clinically meaningful
            prior(cauchy(0,0.5), class = sd)) 

meta.ldl.brm <- brm(SMD|se(Pooled.LDL.SD) ~ 1 + (1|Study),
             data = meta.data,
             prior = priors,
             iter = 4000,
             control=list(adapt_delta=0.9)) #to fix warning There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

```
## 
## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 4.2e-05 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.42 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 4000 [  0%]  (Warmup)
## Chain 1: Iteration:  400 / 4000 [ 10%]  (Warmup)
## Chain 1: Iteration:  800 / 4000 [ 20%]  (Warmup)
## Chain 1: Iteration: 1200 / 4000 [ 30%]  (Warmup)
## Chain 1: Iteration: 1600 / 4000 [ 40%]  (Warmup)
## Chain 1: Iteration: 2000 / 4000 [ 50%]  (Warmup)
## Chain 1: Iteration: 2001 / 4000 [ 50%]  (Sampling)
## Chain 1: Iteration: 2400 / 4000 [ 60%]  (Sampling)
## Chain 1: Iteration: 2800 / 4000 [ 70%]  (Sampling)
## Chain 1: Iteration: 3200 / 4000 [ 80%]  (Sampling)
## Chain 1: Iteration: 3600 / 4000 [ 90%]  (Sampling)
## Chain 1: Iteration: 4000 / 4000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 0.118 seconds (Warm-up)
## Chain 1:                0.1 seconds (Sampling)
## Chain 1:                0.218 seconds (Total)
## Chain 1: 
## 
## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
## Chain 2: 
## Chain 2: Gradient evaluation took 9e-06 seconds
## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
## Chain 2: Adjust your expectations accordingly!
## Chain 2: 
## Chain 2: 
## Chain 2: Iteration:    1 / 4000 [  0%]  (Warmup)
## Chain 2: Iteration:  400 / 4000 [ 10%]  (Warmup)
## Chain 2: Iteration:  800 / 4000 [ 20%]  (Warmup)
## Chain 2: Iteration: 1200 / 4000 [ 30%]  (Warmup)
## Chain 2: Iteration: 1600 / 4000 [ 40%]  (Warmup)
## Chain 2: Iteration: 2000 / 4000 [ 50%]  (Warmup)
## Chain 2: Iteration: 2001 / 4000 [ 50%]  (Sampling)
## Chain 2: Iteration: 2400 / 4000 [ 60%]  (Sampling)
## Chain 2: Iteration: 2800 / 4000 [ 70%]  (Sampling)
## Chain 2: Iteration: 3200 / 4000 [ 80%]  (Sampling)
## Chain 2: Iteration: 3600 / 4000 [ 90%]  (Sampling)
## Chain 2: Iteration: 4000 / 4000 [100%]  (Sampling)
## Chain 2: 
## Chain 2:  Elapsed Time: 0.126 seconds (Warm-up)
## Chain 2:                0.143 seconds (Sampling)
## Chain 2:                0.269 seconds (Total)
## Chain 2: 
## 
## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
## Chain 3: 
## Chain 3: Gradient evaluation took 1e-05 seconds
## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.1 seconds.
## Chain 3: Adjust your expectations accordingly!
## Chain 3: 
## Chain 3: 
## Chain 3: Iteration:    1 / 4000 [  0%]  (Warmup)
## Chain 3: Iteration:  400 / 4000 [ 10%]  (Warmup)
## Chain 3: Iteration:  800 / 4000 [ 20%]  (Warmup)
## Chain 3: Iteration: 1200 / 4000 [ 30%]  (Warmup)
## Chain 3: Iteration: 1600 / 4000 [ 40%]  (Warmup)
## Chain 3: Iteration: 2000 / 4000 [ 50%]  (Warmup)
## Chain 3: Iteration: 2001 / 4000 [ 50%]  (Sampling)
## Chain 3: Iteration: 2400 / 4000 [ 60%]  (Sampling)
## Chain 3: Iteration: 2800 / 4000 [ 70%]  (Sampling)
## Chain 3: Iteration: 3200 / 4000 [ 80%]  (Sampling)
## Chain 3: Iteration: 3600 / 4000 [ 90%]  (Sampling)
## Chain 3: Iteration: 4000 / 4000 [100%]  (Sampling)
## Chain 3: 
## Chain 3:  Elapsed Time: 0.12 seconds (Warm-up)
## Chain 3:                0.109 seconds (Sampling)
## Chain 3:                0.229 seconds (Total)
## Chain 3: 
## 
## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
## Chain 4: 
## Chain 4: Gradient evaluation took 8e-06 seconds
## Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
## Chain 4: Adjust your expectations accordingly!
## Chain 4: 
## Chain 4: 
## Chain 4: Iteration:    1 / 4000 [  0%]  (Warmup)
## Chain 4: Iteration:  400 / 4000 [ 10%]  (Warmup)
## Chain 4: Iteration:  800 / 4000 [ 20%]  (Warmup)
## Chain 4: Iteration: 1200 / 4000 [ 30%]  (Warmup)
## Chain 4: Iteration: 1600 / 4000 [ 40%]  (Warmup)
## Chain 4: Iteration: 2000 / 4000 [ 50%]  (Warmup)
## Chain 4: Iteration: 2001 / 4000 [ 50%]  (Sampling)
## Chain 4: Iteration: 2400 / 4000 [ 60%]  (Sampling)
## Chain 4: Iteration: 2800 / 4000 [ 70%]  (Sampling)
## Chain 4: Iteration: 3200 / 4000 [ 80%]  (Sampling)
## Chain 4: Iteration: 3600 / 4000 [ 90%]  (Sampling)
## Chain 4: Iteration: 4000 / 4000 [100%]  (Sampling)
## Chain 4: 
## Chain 4:  Elapsed Time: 0.122 seconds (Warm-up)
## Chain 4:                0.111 seconds (Sampling)
## Chain 4:                0.233 seconds (Total)
## Chain 4:
```

``` r
plot(meta.ldl.brm, ask=F)
```

![](figures/ldlc-bayesian-meta-1.png)<!-- -->

``` r
summary(meta.ldl.brm)
```

```
##  Family: gaussian 
##   Links: mu = identity; sigma = identity 
## Formula: SMD | se(Pooled.LDL.SD) ~ 1 + (1 | Study) 
##    Data: meta.data (Number of observations: 11) 
##   Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
##          total post-warmup draws = 8000
## 
## Multilevel Hyperparameters:
## ~Study (Number of levels: 11) 
##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sd(Intercept)     1.15      2.16     0.02     7.18 1.00    12545     4464
## 
## Regression Coefficients:
##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## Intercept    18.05      4.59     8.89    27.10 1.00    13910     5788
## 
## Further Distributional Parameters:
##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sigma     0.00      0.00     0.00     0.00   NA       NA       NA
## 
## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).
```

``` r
pp_check(meta.ldl.brm)
```

![](figures/ldlc-bayesian-meta-2.png)<!-- -->

``` r
ranef(meta.ldl.brm)
```

```
## $Study
## , , Intercept
## 
##               Estimate Est.Error  Q2.5 Q97.5
## Al-Sarraj2010  -0.0341      2.31 -3.44  3.47
## Goday2016      -0.0581      2.39 -4.17  3.46
## Hyde2019       -0.0759      2.34 -3.91  3.39
## Saslow2014     -0.0984      2.59 -4.23  3.44
## Sharman2002    -0.0150      2.35 -3.79  3.49
## Sharman2004    -0.0460      2.50 -3.93  3.50
## Sun2019        -0.0138      2.43 -3.78  3.80
## Urbain2017      0.0204      2.45 -3.78  3.66
## Volek2003-4w   -0.0185      2.26 -3.51  3.66
## Volek2009      -0.0167      2.33 -3.90  3.60
## Volek2013      -0.0377      2.40 -3.71  3.30
```

``` r
post.samples <- as_draws_df(meta.ldl.brm) %>%
  select(b_Intercept,sd_Study__Intercept) %>%
  rename(tau=sd_Study__Intercept,
         smd=b_Intercept)

# density plots of posterior distributions
ggplot(aes(x = smd), data = post.samples) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samples$smd)) +
  labs(x = expression(italic(SMD)),
       y = element_blank()) +
  theme_minimal()
```

![](figures/ldlc-bayesian-meta-3.png)<!-- -->

``` r
ggplot(aes(x = tau), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(post.samples$tau)) +        # add point at mean
    labs(x = expression(tau),
       y = element_blank()) +
  theme_minimal()
```

![](figures/ldlc-bayesian-meta-4.png)<!-- -->

``` r
#calculating exact
smd.ecdf <- ecdf(post.samples$smd)
smd.ecdf(0) #probability that effect is greater than zero
```

```
## [1] 0.00025
```

``` r
library(tidybayes)
library(dplyr)
library(ggplot2)
library(ggridges)
library(glue)
library(stringr)
library(forcats)

study.draws <- spread_draws(meta.ldl.brm, r_Study[Study,], b_Intercept) %>% 
  mutate(b_Intercept = r_Study + b_Intercept)

pooled.effect.draws <- spread_draws(meta.ldl.brm, b_Intercept) %>% 
  mutate(Study = "Pooled Effect")

forest.data.bayes <- bind_rows(study.draws, 
                         pooled.effect.draws) %>% 
   ungroup() %>%
   mutate(Study = str_replace_all(Study, "[.]", " ")) %>% 
   mutate(Study = reorder(Study, b_Intercept))

forest.data.summary <- group_by(forest.data.bayes, Study) %>% 
  mean_qi(b_Intercept)

ggplot(aes(b_Intercept, 
           relevel(Study, "Pooled Effect", 
                   after = Inf)), 
       data = forest.data.bayes) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(meta.ldl.brm)[1, 1], 
             color = "grey", size = 1) +
  geom_vline(xintercept = fixef(meta.ldl.brm)[1, 3:4], 
             color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black", 
             size = 1) +
  
  # Add densities
  geom_density_ridges(fill = "blue", 
                      rel_min_height = 0.01, 
                      col = NA, scale = 1,
                      alpha = 0.8) +
  geom_pointintervalh(data = forest.data.summary, 
                      size = 1) +
  
  # Add text and labels
  geom_text(data = mutate_if(forest.data.summary, 
                             is.numeric, round, 2),
    aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"), 
        x = Inf), hjust = "inward") +
  labs(x = "Standardized Mean Difference", # summary measure
       y = element_blank()) +
  lims(x=c(-5,5)) +
  theme_minimal()
```

![](figures/ldlc-bayesian-meta-5.png)<!-- -->
=======

### Body Weight


``` r
weight.meta.data <-
  eval.data %>%
  mutate(`Weight Change SD` = case_when(is.na(`Weight Change SD`)~`Weight Change CI`/1.96,
                                      TRUE~`Weight Change SD`)) %>% 
  #dplyr::select(`Study`,`Weight Change SD`, `Weight Change CI`, `Weight Change`,`Baseline Weight SD`)
  filter(!is.na(`Weight Change SD`)) %>%
  #filter(is.na(`Change in LDL-C`)) %>%
  mutate(MD.Wt = `Weight Change`) %>%
  mutate(BMD=`BMI Change`)

weight.change.meta <- metagen(TE = MD.Wt,
                 seTE = `Weight Change SD`,
                 n.e=n,
                 n.c=n,
                 studlab = Study,
                 data = weight.meta.data,
                 sm = "MD",
                 comb.fixed = TRUE,
                 comb.random = FALSE,
                 #method.tau = "REML",
                 hakn = TRUE,
                 title = "Weight Changes in Ketogenic Diet Studies")

forest(weight.change.meta, 
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

![](figures/weight-meta--1.png)<!-- -->

We evaluated 16 studies for this meta-analysis. Using the meta-analysis method, we found fasting blood LDL-C levels were increased 7.667 mg/dL (95% CI: -3.932 to 19.267) after the ketogenic diet intervention compared to pre-intervention levels, with a significant p-value of 0.195. Across these studies, the I<sup>2</sup> is 0, the p-value for Q is 0.997. This is a highly consistent I^2. 
>>>>>>> Stashed changes

# Average Change in LDL-C


``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=reorder(Study, -`Change in LDL-C`))) +
  geom_bar(stat='identity') +
  labs(x="Study")
```

![](figures/ldl-change-1.png)<!-- -->

``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=reorder(Study, -`Change in LDL-C`))) +
  geom_bar(stat='identity', aes(fill=`Percent Male`)) +
  labs(x="Study")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](figures/ldl-change-2.png)<!-- -->

``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=reorder(Study, -`Change in LDL-C`))) +
  geom_bar(stat='identity', aes(fill=`Sex.Group`)) +
  labs(x="Study") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](figures/ldl-change-3.png)<!-- -->

## Relative to Weight


``` r
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

``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline BMI`)) +
  geom_point() +
  geom_smooth(span=1) +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline BMI")
```

![](figures/ldl-change-vs-weight-2.png)<!-- -->

``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline Weight`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline Weight")
```

![](figures/ldl-change-vs-weight-3.png)<!-- -->

``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`Baseline BMI`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="Baseline BMI")
```

![](figures/ldl-change-vs-weight-4.png)<!-- -->

``` r
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

``` r
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

``` r
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

``` r
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

``` r
library(broom)
bind_rows(`Baseline`=shapiro.test(eval.data$`Baseline Weight`)$p.value,
          `Change`=shapiro.test(eval.data$`Change in LDL-C`)$p.value) %>%
  kable(caption="Shapiro Tests for Correlates")
```



Table: Shapiro Tests for Correlates

| Baseline| Change|
|--------:|------:|
|    0.203|   0.01|

``` r
with(eval.data, cor.test(`Change in LDL-C`,`Baseline Weight`, method="spearman")) %>% tidy %>% kable(caption="Correlation between change in LDL-C and baseline weight")
```



Table: Correlation between change in LDL-C and baseline weight

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|   -0.569|      1789|   0.011|Spearman's rank correlation rho |two.sided   |

``` r
lm(`Change in LDL-C`~`Baseline Weight`+`Percent Male`, data=eval.data) %>% tidy %>% kable(caption="Linear model between change in LDL-C and baseline weight, including gender")
```



Table: Linear model between change in LDL-C and baseline weight, including gender

|term              | estimate| std.error| statistic| p.value|
|:-----------------|--------:|---------:|---------:|-------:|
|(Intercept)       |   42.589|     14.19|     3.001|   0.008|
|`Baseline Weight` |   -0.392|      0.18|    -2.179|   0.045|
|`Percent Male`    |    0.472|     10.91|     0.043|   0.966|

``` r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight`*`Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term                        | df|  sumsq|  meansq| statistic| p.value|
|:---------------------------|--:|------:|-------:|---------:|-------:|
|`Baseline Weight`           |  1| 1081.2| 1081.24|     5.811|   0.031|
|Sex.Group                   |  2|  146.0|   72.98|     0.392|   0.683|
|`Baseline Weight`:Sex.Group |  2|   16.9|    8.47|     0.046|   0.956|
|Residuals                   | 13| 2419.0|  186.08|        NA|      NA|

``` r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  1081|   1081|     6.658|   0.021|
|Sex.Group         |  2|   146|     73|     0.449|   0.646|
|Residuals         | 15|  2436|    162|        NA|      NA|

``` r
ldl.weight.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline Weight`, data = eval.data)
ldl.weight.baseline.lm %>% tidy %>% kable
```



|term              | estimate| std.error| statistic| p.value|
|:-----------------|--------:|---------:|---------:|-------:|
|(Intercept)       |   42.386|    12.999|      3.26|   0.005|
|`Baseline Weight` |   -0.388|     0.145|     -2.67|   0.016|

``` r
ldl.weightsex.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |   43.559|     14.27|     3.053|   0.008|
|`Baseline Weight`    |   -0.412|      0.19|    -2.165|   0.047|
|Sex.GroupMixed       |    2.873|      8.14|     0.353|   0.729|
|Sex.GroupMostly Male |   -6.274|     11.40|    -0.551|   0.590|

``` r
ldl.weightsex.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline Weight` + `Sex.Group`, data = eval.data)
ldl.weightsex.baseline.aov %>% tidy %>% kable
```



|term              | df| sumsq| meansq| statistic| p.value|
|:-----------------|--:|-----:|------:|---------:|-------:|
|`Baseline Weight` |  1|  1081|   1081|     6.658|   0.021|
|Sex.Group         |  2|   146|     73|     0.449|   0.646|
|Residuals         | 15|  2436|    162|        NA|      NA|

``` r
ldl.bmi.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, data = eval.data)
ldl.bmi.baseline.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    42.09|    15.540|      2.71|   0.018|
|`Baseline BMI` |    -1.12|     0.493|     -2.28|   0.040|

``` r
ldl.bmi.baseline.sex.aov <- aov(`Change in LDL-C` ~ `Baseline BMI`*`Sex.Group`, data = eval.data)
ldl.bmi.baseline.sex.aov %>% tidy %>% kable
```



|term                     | df|   sumsq| meansq| statistic| p.value|
|:------------------------|--:|-------:|------:|---------:|-------:|
|`Baseline BMI`           |  1|  709.76| 709.76|     4.491|   0.060|
|Sex.Group                |  2|  193.84|  96.92|     0.613|   0.561|
|`Baseline BMI`:Sex.Group |  1|    5.92|   5.92|     0.037|   0.850|
|Residuals                | 10| 1580.29| 158.03|        NA|      NA|

``` r
ldl.bmi.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline BMI` + `Sex.Group`, data = eval.data)
ldl.bmi.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Baseline BMI` |  1|   710|  709.8|     4.922|   0.048|
|Sex.Group      |  2|   194|   96.9|     0.672|   0.530|
|Residuals      | 11|  1586|  144.2|        NA|      NA|

``` r
ldl.bmi.baseline.mixed.lm<- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mixed"))
ldl.bmi.baseline.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mixed groups")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mixed groups

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   32.463|    10.667|      3.04|   0.019|
|`Baseline BMI` |   -0.839|     0.311|     -2.69|   0.031|

``` r
ldl.bmi.baseline.female.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mostly Female"))
ldl.bmi.baseline.female.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mostly Female groups")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mostly Female groups

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   29.600|     57.93|     0.511|   0.645|
|`Baseline BMI` |   -0.528|      2.31|    -0.229|   0.834|

``` r
ldl.bmi.baseline.male.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, Sex.Group == "Mostly Male"))
ldl.bmi.baseline.male.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Mostly Male")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Mostly Male

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    -7.74|       NaN|       NaN|     NaN|
|`Baseline BMI` |       NA|        NA|        NA|      NA|

``` r
ldl.bmi.baseline.sex.aov <- aov(`Change in LDL-C` ~ `Baseline BMI`*`Normal weight`, data = eval.data)
ldl.bmi.baseline.sex.aov %>% tidy %>% kable
```



|term                           | df|  sumsq| meansq| statistic| p.value|
|:------------------------------|--:|------:|------:|---------:|-------:|
|`Baseline BMI`                 |  1|  709.8|  709.8|     1.328|   0.368|
|`Normal weight`                |  9|  663.8|   73.8|     0.138|   0.987|
|`Baseline BMI`:`Normal weight` |  2|   47.5|   23.8|     0.044|   0.957|
|Residuals                      |  2| 1068.7|  534.3|        NA|      NA|

``` r
ldl.bmi.baseline.female.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, `Normal weight` == "Yes"))
ldl.bmi.baseline.female.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Normal weight studies")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Normal weight studies

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    67.09|       428|     0.157|   0.901|
|`Baseline BMI` |    -2.45|        19|    -0.129|   0.918|

``` r
ldl.bmi.baseline.male.lm <- lm(`Change in LDL-C` ~ `Baseline BMI`, filter(eval.data, `Normal weight` == "No"))
ldl.bmi.baseline.male.lm %>% tidy %>% kable (caption="Correlation betwteen baseline BMI and delta-LDL of Overweight studies")
```



Table: Correlation betwteen baseline BMI and delta-LDL of Overweight studies

|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   -46.76|    104.15|    -0.449|   0.731|
|`Baseline BMI` |     1.36|      2.96|     0.458|   0.727|

Lower baseline BMI was associated with an increased change in LDL-C after consumption of a ketogenic diet (r<sup>2</sup> = 0.285, p-value = 0.04). The association with increased LDL-C was consistent with baseline weight, where a lower baseline weight was associated with an increased change in LDL-C (r<sup>2</sup> = 0.295, p-value = 0.016). 

## Relative to Weight Loss


``` r
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

``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`BMI Change`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="BMI Change")
```

![](figures/ldl-change-vs-weight-change-2.png)<!-- -->

``` r
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

``` r
eval.data %>%
  ggplot(aes(y=`Change in LDL-C`,
             x=`PCT.BMI.Change`)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label=Study)) +
  labs(x="% BMI Change")
```

![](figures/ldl-change-vs-weight-change-4.png)<!-- -->

``` r
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

``` r
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

``` r
library(broom)
bind_rows(`Weight Change`=shapiro.test(eval.data$Pct.Wt.Change)$p.value,
          `LDL Change`=shapiro.test(eval.data$`Change in LDL-C`)$p.value) %>%
  kable(caption="Shapiro Tests for Correlates")
```



Table: Shapiro Tests for Correlates

| Weight Change| LDL Change|
|-------------:|----------:|
|         0.076|       0.01|

``` r
with(eval.data, cor.test(`Change in LDL-C`,Pct.Wt.Change, method="spearman")) %>% tidy %>% kable(caption="Correlation between change in LDL-C and weight change")
```



Table: Correlation between change in LDL-C and weight change

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|    0.444|       539|   0.065|Spearman's rank correlation rho |two.sided   |

``` r
ldl.weightsexinter.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change`*`Sex.Group`, data = eval.data)
ldl.weightsexinter.change.aov %>% tidy %>% kable
```



|term                    | df|  sumsq| meansq| statistic| p.value|
|:-----------------------|--:|------:|------:|---------:|-------:|
|Pct.Wt.Change           |  1|  632.9|  632.9|     3.930|   0.071|
|Sex.Group               |  2|  609.4|  304.7|     1.892|   0.193|
|Pct.Wt.Change:Sex.Group |  2|   92.7|   46.4|     0.288|   0.755|
|Residuals               | 12| 1932.4|  161.0|        NA|      NA|

``` r
ldl.weightsex.change.aov <- aov(`Change in LDL-C` ~ `Pct.Wt.Change` + `Sex.Group`, data = eval.data)
ldl.weightsex.change.aov %>% tidy %>% kable
```



|term          | df| sumsq| meansq| statistic| p.value|
|:-------------|--:|-----:|------:|---------:|-------:|
|Pct.Wt.Change |  1|   633|    633|      4.38|   0.055|
|Sex.Group     |  2|   609|    305|      2.11|   0.159|
|Residuals     | 14|  2025|    145|        NA|      NA|

``` r
ldl.weight.change.lm <- lm(`Change in LDL-C` ~ `Pct.Wt.Change`, data = eval.data)
ldl.weight.change.lm %>% tidy %>% kable
```



|term          | estimate| std.error| statistic| p.value|
|:-------------|--------:|---------:|---------:|-------:|
|(Intercept)   |    19.39|     5.833|      3.32|   0.004|
|Pct.Wt.Change |     1.56|     0.798|      1.96|   0.068|

``` r
ldl.bmisexinter.change.aov <- aov(`Change in LDL-C` ~ `BMI Change`*`Sex.Group`, data = eval.data)
ldl.bmisexinter.change.aov %>% tidy %>% kable
```



|term                   | df| sumsq| meansq| statistic| p.value|
|:----------------------|--:|-----:|------:|---------:|-------:|
|`BMI Change`           |  1|   510|  510.0|     10.05|   0.025|
|Sex.Group              |  1|   590|  589.9|     11.62|   0.019|
|`BMI Change`:Sex.Group |  1|   386|  386.2|      7.61|   0.040|
|Residuals              |  5|   254|   50.8|        NA|      NA|

``` r
ldl.bmisex.change.aov <- aov(`Change in LDL-C` ~ `BMI Change` + `Sex.Group`, data = eval.data)
ldl.bmisex.change.aov %>% tidy %>% kable
```



|term         | df| sumsq| meansq| statistic| p.value|
|:------------|--:|-----:|------:|---------:|-------:|
|`BMI Change` |  1|   510|    510|      4.78|   0.071|
|Sex.Group    |  1|   590|    590|      5.53|   0.057|
|Residuals    |  6|   640|    107|        NA|      NA|

``` r
ldl.bmi.change.lm <- lm(`Change in LDL-C` ~ `BMI Change`, data = eval.data)
ldl.bmi.change.lm %>% tidy %>% kable
```



|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |     22.4|      8.46|      2.64|   0.033|
|`BMI Change` |      5.1|      2.99|      1.70|   0.132|

``` r
ldl.bmi.change.mixed.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mixed"))
ldl.bmi.change.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mixed groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mixed groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |   -0.611|      6.28|    -0.097|   0.927|
|`BMI Change` |   -0.952|      1.88|    -0.506|   0.639|

``` r
#ldl.bmi.change.male.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Male"))
#ldl.bmi.change.male.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Male groups")

<<<<<<< Updated upstream
ldl.bmi.change.female.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Female"))
ldl.bmi.change.female.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Female groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mostly Female groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |     46.5|     11.91|      3.90|   0.160|
|`BMI Change` |     18.5|      9.27|      1.99|   0.296|

``` r
ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
ldl.pctbmi.change.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    20.17|      9.18|      2.20|   0.064|
|PCT.BMI.Change |     1.37|      1.06|      1.28|   0.240|

``` r
ldl.pctbmisexinter.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change`*`Sex.Group`, data = eval.data)
ldl.pctbmisexinter.change.aov %>% tidy %>% kable
```



|term                     | df| sumsq| meansq| statistic| p.value|
|:------------------------|--:|-----:|------:|---------:|-------:|
|PCT.BMI.Change           |  1|   331|  331.4|      4.92|   0.077|
|Sex.Group                |  1|   760|  760.0|     11.28|   0.020|
|PCT.BMI.Change:Sex.Group |  1|   311|  311.4|      4.62|   0.084|
|Residuals                |  5|   337|   67.4|        NA|      NA|

``` r
ldl.pctbmisex.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change` + `Sex.Group`, data = eval.data)
ldl.pctbmisex.change.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|PCT.BMI.Change |  1|   331|    331|      3.07|   0.130|
|Sex.Group      |  1|   760|    760|      7.03|   0.038|
|Residuals      |  6|   648|    108|        NA|      NA|

``` r
ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
ldl.pctbmi.change.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |    20.17|      9.18|      2.20|   0.064|
|PCT.BMI.Change |     1.37|      1.06|      1.28|   0.240|

``` r
ldl.pctbmi.change.mixed.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mixed"))
ldl.pctbmi.change.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mixed groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mixed groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |   -0.611|      6.28|    -0.097|   0.927|
|`BMI Change` |   -0.952|      1.88|    -0.506|   0.639|

``` r
#ldl.pctbmi.change.male.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Male"))
#ldl.pctbmi.change.male.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Male groups")

ldl.pctbmi.change.female.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Female"))
ldl.pctbmi.change.female.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Female groups")
```



Table: Correlation betwteen BMI Change and delta-LDL of Mostly Female groups

|term         | estimate| std.error| statistic| p.value|
|:------------|--------:|---------:|---------:|-------:|
|(Intercept)  |     46.5|     11.91|      3.90|   0.160|
|`BMI Change` |     18.5|      9.27|      1.99|   0.296|
Greater BMI decreases over the study period were associated with a smaller increase in LDL-C after consumption of a ketogenic diet, though this did not reach significance (r<sup>2</sup> = 0.293, p-value = 0.132). The association with the change in LDL-C and decrease in BMI was consistent with weight, with change in weight on LDL-C reaching significance (r<sup>2</sup> = 0.194, p-value = 0.068), where greater decreases in weight were associated with lower increases in LDL-C after consumption of a ketogenic diet. Looking at percent BMI change to account for baseline BMI, greater percent change decreases were associated with a lower increase in LDL-C on a ketogenic diet, though this was not significant (r<sup>2</sup> = 0.19, p-value = 0.24).
=======
# ldl.bmi.change.female.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Female"))
# ldl.bmi.change.female.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Female groups")
# 
# ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
# ldl.pctbmi.change.lm %>% tidy %>% kable
# 
# ldl.pctbmisexinter.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change`*`Sex.Group`, data = eval.data)
# ldl.pctbmisexinter.change.aov %>% tidy %>% kable
# 
# ldl.pctbmisex.change.aov <- aov(`Change in LDL-C` ~ `PCT.BMI.Change` + `Sex.Group`, data = eval.data)
# ldl.pctbmisex.change.aov %>% tidy %>% kable
# 
# ldl.pctbmi.change.lm <- lm(`Change in LDL-C` ~ `PCT.BMI.Change`, data = eval.data)
# ldl.pctbmi.change.lm %>% tidy %>% kable
# 
# ldl.pctbmi.change.mixed.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mixed"))
# ldl.pctbmi.change.mixed.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mixed groups")
# 
# ldl.pctbmi.change.male.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Male"))
# ldl.pctbmi.change.male.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Male groups")
# 
# ldl.pctbmi.change.female.lm<- lm(`Change in LDL-C` ~ `BMI Change`, filter(eval.data, Sex.Group == "Mostly Female"))
# ldl.pctbmi.change.female.lm %>% tidy %>% kable (caption="Correlation betwteen BMI Change and delta-LDL of Mostly Female groups")
```
Greater BMI decreases over the study period were associated with a smaller increase in LDL-C after consumption of a ketogenic diet, though this did not reach significance (r<sup>2</sup> = ` ldl.bmi.change.lm %>% glance %>% pull(r.squared)`, p-value = ` ldl.bmi.change.lm %>% tidy %>% pull(p.value) %>% last`). The association with the change in LDL-C and decrease in BMI was consistent with weight, with change in weight on LDL-C reaching significance (r<sup>2</sup> = ` ldl.weight.change.lm %>% glance %>% pull(r.squared)`, p-value = ` ldl.weight.change.lm %>% tidy %>% pull(p.value) %>% last`), where greater decreases in weight were associated with lower increases in LDL-C after consumption of a ketogenic diet. Looking at percent BMI change to account for baseline BMI, greater percent change decreases were associated with a lower increase in LDL-C on a ketogenic diet, though this was not significant (r<sup>2</sup> = ` ldl.pctbmi.change.lm %>% glance %>% pull(r.squared)`, p-value = ` ldl.pctbmi.change.lm %>% tidy %>% pull(p.value) %>% last`).
>>>>>>> Stashed changes


## Relative to Baseline LDL-C


``` r
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

``` r
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

``` r
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

``` r
library(broom)
bind_rows(`Baseline`=shapiro.test(eval.data$`Baseline LDL`)$p.value,
          `Change`=shapiro.test(eval.data$`Change in LDL-C`)$p.value) %>%
  kable(caption="Shapiro Tests for Correlates")
```



Table: Shapiro Tests for Correlates

| Baseline| Change|
|--------:|------:|
|    0.143|   0.01|

``` r
with(eval.data, cor.test(`Change in LDL-C`,`Baseline LDL`, method="spearman")) %>% tidy %>% kable(caption="Correlation between baseline and change in LDL-C")
```



Table: Correlation between baseline and change in LDL-C

| estimate| statistic| p.value|method                          |alternative |
|--------:|---------:|-------:|:-------------------------------|:-----------|
|    0.011|      1315|   0.962|Spearman's rank correlation rho |two.sided   |

``` r
ldl.baselinesexinter.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` * `Sex.Group`, data = eval.data)
ldl.baselinesexinter.aov %>% tidy %>% kable
```



|term                     | df|   sumsq| meansq| statistic| p.value|
|:------------------------|--:|-------:|------:|---------:|-------:|
|`Baseline LDL`           |  1|    5.56|   5.56|     0.025|   0.876|
|Sex.Group                |  2|  496.12| 248.06|     1.122|   0.353|
|`Baseline LDL`:Sex.Group |  2|  175.33|  87.66|     0.396|   0.680|
|Residuals                | 14| 3095.59| 221.11|        NA|      NA|

``` r
ldl.baselinesex.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baselinesex.aov %>% tidy %>% kable
```



|term           | df|   sumsq| meansq| statistic| p.value|
|:--------------|--:|-------:|------:|---------:|-------:|
|`Baseline LDL` |  1|    5.56|   5.56|     0.027|   0.871|
|Sex.Group      |  2|  496.12| 248.06|     1.213|   0.323|
|Residuals      | 16| 3270.91| 204.43|        NA|      NA|

``` r
ldl.baseline.lm <- lm(`Change in LDL-C` ~ `Baseline LDL`, data = eval.data)
ldl.baseline.lm %>% tidy %>% kable
```



|term           | estimate| std.error| statistic| p.value|
|:--------------|--------:|---------:|---------:|-------:|
|(Intercept)    |   11.129|    19.484|     0.571|   0.575|
|`Baseline LDL` |   -0.028|     0.169|    -0.163|   0.872|

``` r
ldl.baselinesex.lm <- lm(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baselinesex.lm %>% tidy %>% kable
```



|term                 | estimate| std.error| statistic| p.value|
|:--------------------|--------:|---------:|---------:|-------:|
|(Intercept)          |   12.442|    19.334|     0.644|   0.529|
|`Baseline LDL`       |    0.022|     0.171|     0.128|   0.900|
|Sex.GroupMixed       |   -8.759|     7.252|    -1.208|   0.245|
|Sex.GroupMostly Male |  -16.564|    11.809|    -1.403|   0.180|

``` r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Baseline LDL` + `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df|   sumsq| meansq| statistic| p.value|
|:--------------|--:|-------:|------:|---------:|-------:|
|`Baseline LDL` |  1|    5.56|   5.56|     0.027|   0.871|
|Sex.Group      |  2|  496.12| 248.06|     1.213|   0.323|
|Residuals      | 16| 3270.91| 204.43|        NA|      NA|

``` r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Sex.Group`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term      | df| sumsq| meansq| statistic| p.value|
|:---------|--:|-----:|------:|---------:|-------:|
|Sex.Group |  2|   498|    249|      1.29|     0.3|
|Residuals | 17|  3274|    193|        NA|      NA|

``` r
ldl.baseline.aov <- aov(`Change in LDL-C` ~ `Percent Male`, data = eval.data)
ldl.baseline.aov %>% tidy %>% kable
```



|term           | df| sumsq| meansq| statistic| p.value|
|:--------------|--:|-----:|------:|---------:|-------:|
|`Percent Male` |  1|   321|    321|      1.67|   0.212|
|Residuals      | 18|  3452|    192|        NA|      NA|
<<<<<<< Updated upstream

Among individuals, baseline LDL-C was not positively correlated with change in LDL-C after consumption of a ketogenic diet and the relationship was not significant (r<sup>2</sup> = 0.001, p-value = 0.872). 



``` r
library(broom)
lm(formula=`Change in LDL-C`~`Baseline Weight`,
   data=eval.data)%>%
  summary %>%
  tidy %>%
  kable
```



|term              | estimate| std.error| statistic| p.value|
|:-----------------|--------:|---------:|---------:|-------:|
|(Intercept)       |   42.386|    12.999|      3.26|   0.005|
|`Baseline Weight` |   -0.388|     0.145|     -2.67|   0.016|

``` r
lm(formula=`Change in LDL-C`~`Baseline HOMA-IR`,
   data=eval.data)%>%
  summary %>%
  tidy %>%
  kable
```



|term               | estimate| std.error| statistic| p.value|
|:------------------|--------:|---------:|---------:|-------:|
|(Intercept)        |    3.682|      5.32|     0.692|    0.52|
|`Baseline HOMA-IR` |   -0.831|      1.37|    -0.608|    0.57|

# Decision Tree


``` r
eval.data %>% select(where(is.numeric)) %>%
  select(-ends_with('SD'),-ends_with('CI')) -> combined.data.numeric

correlation.matrix <- cor(combined.data.numeric, use="pairwise.complete.obs")
library(corrplot)
corrplot(correlation.matrix, method="color", diag=FALSE, order="alphabet")
```

![](figures/decision-tree-1.png)<!-- -->

``` r
corrplot(correlation.matrix, method="color", diag=FALSE, tl.cex = .7)
```

![](figures/decision-tree-2.png)<!-- -->

``` r
library(rpart)
combined.data.numeric -> tree.data
rpart(`Change in LDL-C`~., data=tree.data) ->tree
library(rattle)
library(rpart.plot)
rpart.plot(tree)
```

![](figures/decision-tree-3.png)<!-- -->


=======

Among individuals, baseline LDL-C was not positively correlated with change in LDL-C after consumption of a ketogenic diet and the relationship was not significant (r<sup>2</sup> = 0.001, p-value = 0.872). 
>>>>>>> Stashed changes

# Session Information


``` r
sessionInfo()
```

```
<<<<<<< Updated upstream
## R version 4.4.1 (2024-06-14)
## Platform: x86_64-apple-darwin20
## Running under: macOS Monterey 12.7.5
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
=======
## R version 4.5.1 (2025-06-13)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.7.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
>>>>>>> Stashed changes
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Detroit
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
<<<<<<< Updated upstream
##  [1] rpart.plot_3.1.2 rattle_5.5.1     bitops_1.0-8     tibble_3.2.1    
##  [5] rpart_4.1.23     corrplot_0.92    broom_1.0.6      ggrepel_0.9.5   
##  [9] forcats_1.0.0    stringr_1.5.1    glue_1.7.0       ggridges_0.5.6  
## [13] tidybayes_3.0.6  brms_2.21.0      Rcpp_1.0.13      meta_7.0-0      
## [17] metadat_1.2-0    readr_2.1.5      ggplot2_3.5.1    dplyr_1.1.4     
## [21] tidyr_1.3.1      knitr_1.48      
## 
## loaded via a namespace (and not attached):
##  [1] gridExtra_2.3        inline_0.3.19        rlang_1.1.4         
##  [4] magrittr_2.0.3       matrixStats_1.3.0    compiler_4.4.1      
##  [7] mgcv_1.9-1           loo_2.8.0            callr_3.7.6         
## [10] vctrs_0.6.5          reshape2_1.4.4       pkgconfig_2.0.3     
## [13] arrayhelpers_1.1-0   crayon_1.5.3         fastmap_1.2.0       
## [16] backports_1.5.0      labeling_0.4.3       utf8_1.2.4          
## [19] rmarkdown_2.27       tzdb_0.4.0           ps_1.7.7            
## [22] nloptr_2.1.1         purrr_1.0.2          bit_4.0.5           
## [25] xfun_0.46            cachem_1.1.0         jsonlite_1.8.8      
## [28] highr_0.11           parallel_4.4.1       R6_2.5.1            
## [31] bslib_0.8.0          stringi_1.8.4        StanHeaders_2.32.10 
## [34] boot_1.3-30          jquerylib_0.1.4      numDeriv_2016.8-1.1 
## [37] rstan_2.32.6         bayesplot_1.11.1     Matrix_1.7-0        
## [40] splines_4.4.1        tidyselect_1.2.1     abind_1.4-5         
## [43] yaml_2.3.10          codetools_0.2-20     metafor_4.6-0       
## [46] processx_3.8.4       pkgbuild_1.4.4       lattice_0.22-6      
## [49] plyr_1.8.9           withr_3.0.0          bridgesampling_1.1-2
## [52] posterior_1.6.0      coda_0.19-4.1        evaluate_0.24.0     
## [55] CompQuadForm_1.4.3   RcppParallel_5.1.8   ggdist_3.3.2        
## [58] xml2_1.3.6           pillar_1.9.0         tensorA_0.36.2.1    
## [61] checkmate_2.3.2      stats4_4.4.1         distributional_0.4.0
## [64] generics_0.1.3       vroom_1.6.5          mathjaxr_1.6-0      
## [67] hms_1.1.3            rstantools_2.4.0     munsell_0.5.1       
## [70] scales_1.3.0         minqa_1.2.7          tools_4.4.1         
## [73] lme4_1.1-35.5        mvtnorm_1.2-5        grid_4.4.1          
## [76] QuickJSR_1.3.1       colorspace_2.1-1     nlme_3.1-164        
## [79] cli_3.6.3            fansi_1.0.6          svUnit_1.0.6        
## [82] Brobdingnag_1.2-9    gtable_0.3.5         sass_0.4.9          
## [85] digest_0.6.36        farver_2.1.2         htmltools_0.5.8.1   
## [88] lifecycle_1.0.4      bit64_4.0.5          MASS_7.3-60.2
=======
## [1] broom_1.0.10  ggrepel_0.9.6 meta_8.2-1    metadat_1.4-0 ggplot2_4.0.0
## [6] readr_2.1.5   dplyr_1.1.4   tidyr_1.3.1   knitr_1.50   
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.6        xfun_0.53           bslib_0.9.0        
##  [4] CompQuadForm_1.4.4  lattice_0.22-7      mathjaxr_1.8-0     
##  [7] tzdb_0.5.0          numDeriv_2016.8-1.1 vctrs_0.6.5        
## [10] tools_4.5.1         Rdpack_2.6.4        generics_0.1.4     
## [13] parallel_4.5.1      tibble_3.3.0        pkgconfig_2.0.3    
## [16] Matrix_1.7-4        RColorBrewer_1.1-3  S7_0.2.0           
## [19] lifecycle_1.0.4     compiler_4.5.1      farver_2.1.2       
## [22] stringr_1.5.2       htmltools_0.5.8.1   sass_0.4.10        
## [25] yaml_2.3.10         pillar_1.11.1       nloptr_2.2.1       
## [28] crayon_1.5.3        jquerylib_0.1.4     MASS_7.3-65        
## [31] cachem_1.1.0        reformulas_0.4.1    boot_1.3-32        
## [34] nlme_3.1-168        tidyselect_1.2.1    digest_0.6.37      
## [37] stringi_1.8.7       purrr_1.1.0         labeling_0.4.3     
## [40] splines_4.5.1       fastmap_1.2.0       grid_4.5.1         
## [43] cli_3.6.5           metafor_4.8-0       magrittr_2.0.4     
## [46] withr_3.0.2         backports_1.5.0     scales_1.4.0       
## [49] bit64_4.6.0-1       rmarkdown_2.30      bit_4.6.0          
## [52] lme4_1.1-37         hms_1.1.4           evaluate_1.0.5     
## [55] rbibutils_2.3       mgcv_1.9-3          rlang_1.1.6        
## [58] Rcpp_1.1.0          glue_1.8.0          xml2_1.4.0         
## [61] rstudioapi_0.17.1   vroom_1.6.6         minqa_1.2.8        
## [64] jsonlite_2.0.0      R6_2.6.1
>>>>>>> Stashed changes
```
