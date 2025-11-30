---
title: "MR Analysis of SNPs for Calcium and their Effect on LDL Cholesterol"
author: "Dave Bridges"
date: "October 14, 2025"
format:
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
  pdf: default
knitr:
  opts_chunk:
    fig.path: "figures/"
    dev: ["png", "pdf"]  # Remove !expr, just use array syntax
    fig.keep: "all"
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


## Purpose

To validate SNPs for calcium GWAS using those identified using UK Biobank.  This script can be found in /Users/davebrid/Documents/GitHub/PrecisionNutrition/Human Genetics and was most recently run on Sun Nov 30 10:39:36 2025

## Data Entry


::: {.cell}

```{.r .cell-code}
instruments.calcium.file <- 'Calcium Instruments from UKBB.csv'
gwas.ldlc.file <- 'PheWeb Summary Statistics/phenocode-LDL.tsv.gz'
samplesize.outcome.ldlc <- 46100


# loaded and renamed columns
instruments.calcium <- read_csv(instruments.calcium.file) |>
  rename(
    SNP                       = SP2,
    beta.exposure             = BETA,
    se.exposure               = SE,
    effect_allele.exposure    = EA,
    other_allele.exposure     = OA,
    pval.exposure             = P,
    eaf.exposure              = ALT_FREQS,
    samplesize.exposure       = N_exposure
  ) |>
  mutate(id.exposure="Calcium (UK Biobank)",
         exposure="Calcium (UK Biobank)")


gwas.ldlc <- read_tsv(gwas.ldlc.file) |>
  mutate(ID=paste(chrom, pos, ref,alt, sep=":")) |>
  rename(
    SNP                        = ID,            # or ID if that’s the matching ID
    beta.outcome               = beta,
    se.outcome                 = sebeta,
    effect_allele.outcome      = alt,   # whichever is effect allele
    other_allele.outcome       = ref,   # whichever is other allele
    pval.outcome               = pval,
    eaf.outcome                = maf,
  ) |>
  mutate(id.outcome = "LDL Cholesterol (MGI-BioVU LabWAS)",
         outcome = "LDL Cholesterol (MGI-BioVU LabWAS)",
         samplesize.outcome = samplesize.outcome.ldlc)  # sample size for MGI/BioVU for calcium)
```
:::


This presumes the sample sizes was 46100 from Table 1 of https://doi.org/10.1371/journal.pgen.1009077.

Loaded in the instruments for calcium from UK Biobank from the datafile Calcium Instruments from UKBB.csv and the GWAS summary statistics for calcium from the datafile PheWeb Summary Statistics/phenocode-LDL.tsv.gz.


::: {.cell}

```{.r .cell-code}
library(TwoSampleMR)

data <- harmonise_data(instruments.calcium, gwas.ldlc, action = 2)

table(data$mr_keep) |>
  kable(caption="Number of SNPs kept for MR analysis")
```

::: {.cell-output-display}


Table: Number of SNPs kept for MR analysis

|Var1  | Freq|
|:-----|----:|
|FALSE |    2|
|TRUE  |  275|


:::

```{.r .cell-code}
table(data$palindromic)  |>
  kable(caption="Number of palindromic SNPs")
```

::: {.cell-output-display}


Table: Number of palindromic SNPs

|Var1  | Freq|
|:-----|----:|
|FALSE |  255|
|TRUE  |   22|


:::

```{.r .cell-code}
data <- data %>%
  mutate(
    allele_match = (toupper(effect_allele.exposure) == toupper(effect_allele.outcome)) &
                  (toupper(other_allele.exposure) == toupper(other_allele.outcome)),
    allele_swapped = (toupper(effect_allele.exposure) == toupper(other_allele.outcome)) &
                    (toupper(other_allele.exposure) == toupper(effect_allele.outcome))
  )

# 2) EAF concordance checks (detect possible strand/orientation issues)
# requires eaf.exposure and eaf.outcome present
if(all(c("eaf.exposure","eaf.outcome") %in% names(data))){
  data <- data %>%
    mutate(
      eaf_diff = abs(eaf.exposure - eaf.outcome),
      eaf_flip_diff = abs(eaf.exposure - (1 - eaf.outcome)),
      suspicious_eaf = (eaf_diff > 0.2 & eaf_flip_diff > 0.2)  # very different frequencies
    )
  summary(data$eaf_diff)
  summary(data$eaf_flip_diff)
  cat("Num suspicious EAFs:", sum(data$suspicious_eaf, na.rm=TRUE), "\n")
} else {
  cat("No EAF columns present for both datasets; consider adding reference panel EAFs.\n")
}
```

::: {.cell-output .cell-output-stdout}

```
Num suspicious EAFs: 0 
```


:::

```{.r .cell-code}
# 3) List discordant SNPs
data <- data %>%
  mutate(sign_match = sign(beta.exposure) == sign(beta.outcome))

discordant <- data %>% filter(!sign_match) %>%
  select(SNP, beta.exposure, se.exposure, beta.outcome, se.outcome,
         effect_allele.exposure, other_allele.exposure,
         effect_allele.outcome, other_allele.outcome,
         palindromic, ambiguous, eaf.exposure, eaf.outcome)

kable(discordant |>
        arrange(beta.exposure) |>
        select(SNP,beta.exposure,se.exposure,beta.outcome,se.outcome),
      caption="Discordant SNPs where the beta coefficients directionally differ between exposure and outcome")
```

::: {.cell-output-display}


Table: Discordant SNPs where the beta coefficients directionally differ between exposure and outcome

|SNP              | beta.exposure| se.exposure| beta.outcome| se.outcome|
|:----------------|-------------:|-----------:|------------:|----------:|
|11:126218773:C:G |      -0.07050|    0.005199|       0.0720|     0.0190|
|2:234264848:C:G  |      -0.04882|    0.002709|       0.0071|     0.0099|
|9:97804641:A:C   |      -0.04572|    0.004387|       0.0250|     0.0170|
|6:74493432:A:C   |      -0.04406|    0.002516|       0.0150|     0.0092|
|4:40555956:A:G   |      -0.04045|    0.004049|       0.0370|     0.0150|
|13:42487165:A:T  |      -0.03960|    0.003880|       0.0110|     0.0140|
|3:12390484:T:C   |      -0.03952|    0.003841|       0.0250|     0.0140|
|10:100179851:T:C |      -0.03719|    0.004556|       0.0011|     0.0160|
|19:50028163:A:G  |      -0.03666|    0.003253|       0.0057|     0.0120|
|9:80366259:T:C   |      -0.03648|    0.002799|       0.0001|     0.0100|
|12:12205282:A:G  |      -0.03404|    0.004502|       0.0008|     0.0170|
|14:95053849:T:C  |      -0.03404|    0.004198|       0.0028|     0.0150|
|11:2956166:T:C   |      -0.03301|    0.002766|       0.0058|     0.0100|
|15:49278090:A:G  |      -0.03252|    0.004474|       0.0190|     0.0170|
|3:121364173:A:T  |      -0.03141|    0.004461|       0.0078|     0.0170|
|1:43458250:T:C   |      -0.03045|    0.002791|       0.0064|     0.0100|
|2:61368532:T:G   |      -0.03039|    0.002531|       0.0021|     0.0093|
|7:65326821:A:G   |      -0.02975|    0.002522|       0.0210|     0.0092|
|1:51659401:A:C   |      -0.02900|    0.005431|       0.0047|     0.0200|
|17:59450441:T:C  |      -0.02881|    0.003364|       0.0170|     0.0120|
|14:92279983:T:G  |      -0.02856|    0.003357|       0.0210|     0.0120|
|14:24862936:T:C  |      -0.02854|    0.005228|       0.0072|     0.0200|
|10:8108592:T:C   |      -0.02730|    0.003026|       0.0200|     0.0110|
|8:106583124:A:G  |      -0.02619|    0.002810|       0.0130|     0.0100|
|8:8283706:A:G    |      -0.02600|    0.003840|       0.0061|     0.0140|
|17:37387413:A:G  |      -0.02579|    0.002984|       0.0280|     0.0110|
|5:131329591:T:C  |      -0.02572|    0.003741|       0.0130|     0.0140|
|17:67290356:A:G  |      -0.02541|    0.004384|       0.0460|     0.0160|
|1:1095130:T:C    |      -0.02524|    0.002786|       0.0094|     0.0110|
|14:94844843:T:G  |      -0.02446|    0.002868|       0.0120|     0.0110|
|7:77433484:A:C   |      -0.02410|    0.002561|       0.0036|     0.0093|
|7:16042596:T:C   |      -0.02389|    0.003055|       0.0210|     0.0110|
|1:155186729:T:C  |      -0.02375|    0.002510|       0.0190|     0.0092|
|2:114032767:T:C  |      -0.02338|    0.002536|       0.0014|     0.0093|
|1:21680008:A:G   |      -0.02306|    0.004535|       0.0110|     0.0170|
|14:105992183:A:C |      -0.02300|    0.002945|       0.0033|     0.0110|
|19:37395697:A:C  |      -0.02241|    0.004269|       0.0078|     0.0160|
|7:92252203:C:G   |      -0.02231|    0.002848|       0.0094|     0.0100|
|7:107202502:A:C  |      -0.02224|    0.002785|       0.0180|     0.0100|
|10:96035980:T:C  |      -0.02221|    0.004032|       0.0200|     0.0150|
|15:78325229:T:C  |      -0.02204|    0.003103|       0.0007|     0.0120|
|3:121323008:C:G  |      -0.02193|    0.003964|       0.0110|     0.0140|
|10:52881906:T:C  |      -0.02148|    0.003597|       0.0060|     0.0130|
|17:47863368:T:C  |      -0.02146|    0.003021|       0.0460|     0.0110|
|10:22337752:A:G  |      -0.02115|    0.002880|       0.0019|     0.0110|
|20:5532853:T:C   |      -0.02072|    0.003083|       0.0094|     0.0120|
|6:137264165:A:G  |      -0.02070|    0.002519|       0.0079|     0.0092|
|2:111989372:T:G  |      -0.02055|    0.004040|       0.0097|     0.0140|
|14:60617639:A:G  |      -0.02018|    0.003536|       0.0093|     0.0130|
|4:27016468:T:C   |      -0.01980|    0.002844|       0.0130|     0.0100|
|2:133188106:T:G  |      -0.01979|    0.003106|       0.0230|     0.0110|
|1:226917898:A:G  |      -0.01971|    0.002520|       0.0100|     0.0092|
|6:25745243:A:T   |      -0.01892|    0.002773|       0.0270|     0.0100|
|21:35893737:T:C  |      -0.01874|    0.002560|       0.0060|     0.0093|
|13:110501545:A:C |      -0.01869|    0.002598|       0.0054|     0.0094|
|17:1982952:A:G   |      -0.01868|    0.002601|       0.0065|     0.0095|
|6:134509254:T:C  |      -0.01851|    0.002664|       0.0310|     0.0098|
|14:70000399:A:C  |      -0.01821|    0.003351|       0.0014|     0.0120|
|9:136154168:T:C  |      -0.01818|    0.003106|       0.0540|     0.0110|
|10:94436851:T:C  |      -0.01808|    0.002518|       0.0044|     0.0092|
|15:96648537:T:C  |      -0.01808|    0.002716|       0.0150|     0.0098|
|19:11979164:T:C  |      -0.01803|    0.002814|       0.0052|     0.0100|
|9:116128568:A:C  |      -0.01762|    0.003277|       0.0110|     0.0120|
|17:17783748:A:G  |      -0.01743|    0.002608|       0.0120|     0.0094|
|3:56865776:A:G   |      -0.01714|    0.002595|       0.0086|     0.0095|
|6:129827116:T:C  |      -0.01691|    0.002776|       0.0120|     0.0100|
|6:44097472:T:C   |      -0.01666|    0.003138|       0.0079|     0.0110|
|7:150614934:T:C  |      -0.01631|    0.002912|       0.0120|     0.0110|
|17:35860964:T:C  |      -0.01609|    0.002611|       0.0091|     0.0095|
|2:190984523:A:T  |      -0.01569|    0.002790|       0.0008|     0.0100|
|11:65363958:A:G  |      -0.01556|    0.002862|       0.0006|     0.0100|
|2:25492467:A:G   |      -0.01555|    0.002718|       0.0098|     0.0098|
|6:139834012:T:G  |      -0.01443|    0.002601|       0.0180|     0.0095|
|10:96305329:T:C  |      -0.01295|    0.002527|       0.0082|     0.0093|
|11:8252853:A:G   |      -0.01279|    0.002444|       0.0057|     0.0092|
|19:14170423:A:G  |       0.01308|    0.002552|      -0.0051|     0.0093|
|7:25997536:A:G   |       0.01348|    0.002555|      -0.0200|     0.0093|
|1:153606851:A:G  |       0.01397|    0.002528|      -0.0059|     0.0092|
|15:60895636:A:G  |       0.01462|    0.002569|      -0.0110|     0.0094|
|1:116804727:A:C  |       0.01493|    0.002574|      -0.0057|     0.0093|
|22:30800338:T:G  |       0.01512|    0.002527|      -0.0091|     0.0092|
|6:43299058:T:G   |       0.01538|    0.002799|      -0.0290|     0.0100|
|8:98812563:A:G   |       0.01540|    0.002523|      -0.0140|     0.0093|
|10:126477406:A:G |       0.01587|    0.002857|      -0.0036|     0.0110|
|15:60762779:T:C  |       0.01602|    0.002763|      -0.0150|     0.0100|
|2:191318776:T:C  |       0.01605|    0.002540|      -0.0057|     0.0093|
|1:200457488:A:G  |       0.01617|    0.002624|      -0.0035|     0.0096|
|21:37835347:A:T  |       0.01630|    0.002865|      -0.0028|     0.0110|
|2:69573029:T:C   |       0.01637|    0.002546|      -0.0052|     0.0092|
|19:30304483:T:C  |       0.01662|    0.002567|      -0.0014|     0.0093|
|12:4057955:T:C   |       0.01683|    0.002947|      -0.0110|     0.0110|
|2:131089135:T:C  |       0.01689|    0.003111|      -0.0082|     0.0110|
|5:55816081:A:C   |       0.01745|    0.002946|      -0.0190|     0.0110|
|10:64914467:A:G  |       0.01789|    0.002722|      -0.0055|     0.0100|
|22:25000229:A:T  |       0.01793|    0.002635|      -0.0041|     0.0096|
|10:104654211:A:G |       0.01805|    0.003067|      -0.0150|     0.0110|
|5:174050806:T:C  |       0.01836|    0.002641|      -0.0020|     0.0097|
|16:81566182:C:G  |       0.01848|    0.002533|      -0.0018|     0.0093|
|19:3133853:T:C   |       0.01856|    0.002842|      -0.0160|     0.0110|
|8:72414019:T:C   |       0.01856|    0.002832|      -0.0049|     0.0100|
|15:52856934:A:C  |       0.01930|    0.003840|      -0.0041|     0.0140|
|1:220009078:A:C  |       0.01973|    0.003802|      -0.0049|     0.0130|
|2:103151136:A:G  |       0.01981|    0.003588|      -0.0014|     0.0130|
|7:143095378:T:C  |       0.01994|    0.003956|      -0.0052|     0.0140|
|1:217456447:T:G  |       0.02116|    0.002513|      -0.0130|     0.0092|
|11:116924978:A:C |       0.02162|    0.003519|      -0.0120|     0.0130|
|13:20436598:T:G  |       0.02175|    0.003751|      -0.0180|     0.0140|
|3:113275624:A:C  |       0.02178|    0.002569|      -0.0075|     0.0093|
|20:33335657:A:G  |       0.02189|    0.003229|      -0.0062|     0.0120|
|12:4573791:T:G   |       0.02200|    0.003770|      -0.0091|     0.0140|
|4:17809538:A:G   |       0.02218|    0.003766|      -0.0240|     0.0140|
|8:21944962:T:G   |       0.02279|    0.002913|      -0.0044|     0.0110|
|8:101703522:A:T  |       0.02295|    0.003722|      -0.0110|     0.0140|
|9:903352:A:T     |       0.02411|    0.004078|      -0.0076|     0.0150|
|5:72385914:A:G   |       0.02412|    0.002674|      -0.0170|     0.0099|
|1:156754465:T:C  |       0.02422|    0.002874|      -0.0055|     0.0110|
|17:73497415:T:C  |       0.02491|    0.004268|      -0.0072|     0.0160|
|14:21944174:A:G  |       0.02525|    0.004435|      -0.0022|     0.0170|
|12:90171438:T:G  |       0.02538|    0.003466|      -0.0130|     0.0130|
|1:21821757:A:G   |       0.02544|    0.002519|      -0.0017|     0.0092|
|2:54846400:A:G   |       0.02574|    0.003572|      -0.0004|     0.0130|
|6:32152387:A:T   |       0.02656|    0.003157|      -0.0110|     0.0110|
|19:3128748:T:C   |       0.02663|    0.003585|      -0.0041|     0.0130|
|17:7731764:T:C   |       0.02673|    0.005343|      -0.0016|     0.0200|
|14:94796184:A:T  |       0.03391|    0.003120|      -0.0084|     0.0110|
|19:19717056:A:G  |       0.03700|    0.004983|      -0.0730|     0.0180|
|3:186337713:T:C  |       0.04188|    0.002626|      -0.0110|     0.0096|
|17:55301369:A:G  |       0.04350|    0.004778|      -0.0250|     0.0180|
|20:5572957:T:C   |       0.04480|    0.006937|      -0.0150|     0.0270|
|3:122069056:T:C  |       0.05004|    0.008926|      -0.0170|     0.0340|
|15:44244523:T:C  |       0.05713|    0.006127|      -0.0250|     0.0220|
|19:35559474:T:C  |       0.05815|    0.003271|      -0.0043|     0.0120|
|3:120753198:T:C  |       0.06249|    0.010840|      -0.0036|     0.0430|
|12:113527940:A:C |       0.07225|    0.011770|      -0.0170|     0.0430|
|1:52787887:T:C   |       0.08565|    0.008950|      -0.0910|     0.0360|
|3:121993247:A:G  |       0.17270|    0.003739|      -0.0120|     0.0130|


:::

```{.r .cell-code}
library(ggrepel)
ggplot(data, aes(x=beta.exposure, y=beta.outcome, color = sign_match)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome), width=0) +
  geom_text_repel(data = filter(data, !sign_match), aes(label=SNP), hjust=0, vjust=0, size=3) +
  theme_minimal() +
  labs(x="beta.exposure", y="beta.outcome", title="Exposure vs Outcome betas; discordant SNPs labeled")
```

::: {.cell-output-display}
![](figures/calcium-ldlc-harmonizing-snps-1.png){width=672}
:::
:::



::: {.cell}

```{.r .cell-code}
ggplot(data, aes(x=beta.exposure, y=beta.outcome)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome,
                    ymax = beta.outcome + 1.96*se.outcome),
                alpha=0.5) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96*se.exposure,
                    xmax = beta.exposure + 1.96*se.exposure),
                alpha=0.5) +
  geom_smooth(method="lm",se=F) +
  theme_classic(base_size=16) +
  labs(x="Exposure Estimate", 
       y="Outcome Estimate", 
       title="") 
```

::: {.cell-output-display}
![](figures/calcium-ldlc-scatter-1.png){width=672}
:::
:::


There were 136 discordant SNPs between the exposure and outcome datasets.  These are listed above.  We can see that some of these SNPs have very small effect sizes in the outcome dataset, suggesting that the discordance may be due to noise.  These were kept in the analysis

### Steiger Filtering


::: {.cell}

```{.r .cell-code}
data_steiger <- steiger_filtering(data)

table(data_steiger$steiger_direction, useNA="ifany") |>
  kable(caption="Steiger filtering results for calcium-LDL evaluation")
```

::: {.cell-output-display}


Table: Steiger filtering results for calcium-LDL evaluation

| Freq|
|----:|


:::
:::


Harmonization results

- We used 363 SNPs as instruments for calcium from UK Biobank.
- There were 275 SNPs in common between the exposure and outcome datasets.
- Removed 0 SNPs due to allele mismatches
- Identified 22 palindromic SNPs 
- A total of 277 SNPs remained for use after harmonization.  2 SNPS were removed because the palindromic SNP is ambiguous and strand alignment could not be resolved, this variant was automatically dropped from the MR analysis to avoid mis-specified effect directions.
- After Steiger filtering, 275 SNPs were retained for analysis, indicating that all SNPs had stronger associations with the exposure (calcium in UK Biobank) than the outcome (calcium in MGI/BioVU), supporting the assumed causal direction.  0 SNPs were removed by Steiger filtering.


::: {.cell}

```{.r .cell-code}
data.annot <- data_steiger %>%
  mutate(
    R2.exposure = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
    F.exposure = (R2.exposure * (samplesize.exposure - 2)) / (1 - R2.exposure)
  )

calcium.exposure.summary <- data.annot %>%
  summarise(
    num_snps = n(),
    samplesize.exposure = first(samplesize.exposure),
    cumulative_R2 = sum(R2.exposure, na.rm = TRUE),
    mean_F = mean(F.exposure, na.rm = TRUE),
    median_F = median(F.exposure, na.rm = TRUE),
    mean_maf = mean(eaf.exposure, na.rm = TRUE),
    mean_beta = mean(abs(beta.exposure), na.rm = TRUE)
  ) |>
  mutate(overall_F = (cumulative_R2 * (samplesize.exposure - num_snps - 1)) / 
                     ((1 - cumulative_R2) * num_snps))

# For outcome (e.g., cholesterol) SNPs
outcome.summary_metrics <- data.annot %>%
  summarise(
    num_snps = n(),
    mean_beta = mean(abs(beta.outcome), na.rm = TRUE),
    mean_se = mean(se.outcome, na.rm = TRUE),
    mean_maf = mean(eaf.outcome, na.rm = TRUE)
  )

library(knitr)
kable(calcium.exposure.summary, caption="Summary of calcium instruments after harmonisation")
```

::: {.cell-output-display}


Table: Summary of calcium instruments after harmonisation

| num_snps| samplesize.exposure| cumulative_R2|  mean_F| median_F|  mean_maf| mean_beta| overall_F|
|--------:|-------------------:|-------------:|-------:|--------:|---------:|---------:|---------:|
|      277|              385066|     0.0638365| 88.8467| 54.59161| 0.3670863| 0.0259874|  94.72378|


:::

```{.r .cell-code}
write_csv(calcium.exposure.summary, "Instrument Metrics - Calcium - Post-Harmonization (LDL-Cholesterol).csv")
write_csv(outcome.summary_metrics, "Outcome Metrics - Calcium - Post-Harmonization (LDL-Cholesterol).csv")

#write out the instruments used for calcium
data_steiger %>% filter(mr_keep==TRUE) %>% 
  mutate(Exposure = "Calcium") |>
  select(Exposure, CHR, POS, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure, R2, `F`, rsids, nearest_genes) |>
  rename(effect_allele = effect_allele.exposure,
         other_allele = other_allele.exposure,
         beta = beta.exposure,
         se = se.exposure,
         p = pval.exposure,
         eaf = eaf.exposure) |>
  write_csv("Calcium Instruments Post-Harmonization (LDL-Cholesterol).csv")
  
kable(outcome.summary_metrics, caption="Summary of LDL cholesterol effects by calcium instruments after harmonisation")
```

::: {.cell-output-display}


Table: Summary of LDL cholesterol effects by calcium instruments after harmonisation

| num_snps| mean_beta|  mean_se|  mean_maf|
|--------:|---------:|--------:|---------:|
|      277| 0.0121996| 0.012583| 0.2610072|


:::
:::


## MR Analyses


::: {.cell}

```{.r .cell-code}
calcium.ldlc.mr <- mr(data_steiger,
                         method_list = c( "mr_ivw_mre",
                                           "mr_ivw_fe",
                                          'mr_raps',
                                         "mr_egger_regression", 
                                         "mr_weighted_median", 
                                         "mr_weighted_mode"))


#M-PRESSO has to be run separately
library(MRPRESSO)

calcium.ldlc.mr_presso_results <- mr_presso(
  BetaOutcome = "beta.outcome",      # Column name for outcome betas
  BetaExposure = "beta.exposure",    # Column name for exposure betas
  SdOutcome = "se.outcome",          # Column name for outcome SEs
  SdExposure = "se.exposure",        # Column name for exposure SEs
  data = data_steiger,                  # Your dataset
  NbDistribution = 2000,              # Number of distributions (default 1000)
  SignifThreshold = 0.05,             # Significance threshold
  OUTLIERtest = TRUE,                 # Perform outlier test
  DISTORTIONtest = TRUE,              # Perform distortion test
)

library(forcats)
calcium.ldlc.mr_presso_results$`Main MR results` |>
  select(`MR Analysis`, `Causal Estimate`, Sd, `P-value`) |>
  rename(method = `MR Analysis`,
         b = `Causal Estimate`,
         se = Sd,
         pval = `P-value`) |>
  mutate(method = fct_recode(method,
                             "MR-PRESSO (Raw)"="Raw",
                             "MR-PRESSO (Outlier-corrected)"="Outlier-corrected")) -> calcium.ldlc.mr_presso_df

calcium.ldlc.mr.mrpresso <- bind_rows(as_tibble(calcium.ldlc.mr), as_tibble(calcium.ldlc.mr_presso_df))|>
  fill(id.exposure, id.outcome, exposure, outcome,nsnp,.direction="down")
```
:::


The primary result, using the inverse variance weighted (multiplicative random effects) method shows a -0.0010879 $\pm$ 0.0383213 SD increase in LDL cholesterol (MGI-BioVU LabWAS) per 1 SD increase in calcium (UK Biobank).  This is statistically significant with a p-value of 0.9773529.  All four MR methods (IVW, weighted mode, weighted median, MR-Egger) gave consistent, **non-significant** causal estimates, supporting the hypothesis of no causal effect of calcium on LDL-C.  

### MR-Egger Intercept


::: {.cell}

```{.r .cell-code}
egger_intercept <- mr_pleiotropy_test(data_steiger)
egger_intercept|>
  select(-starts_with('id')) |> 
  kable(caption="MR Pleiotropy Results for Calcium-Cholesterol Analysis")
```

::: {.cell-output-display}


Table: MR Pleiotropy Results for Calcium-Cholesterol Analysis

|outcome                            |exposure             | egger_intercept|        se|      pval|
|:----------------------------------|:--------------------|---------------:|---------:|---------:|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |       0.0019154| 0.0020949| 0.3613545|


:::
:::

The MR-Egger intercept is  with a p-value of 0.3613545, indicating no evidence of directional pleiotropy.  The intercept magnitude is near zero, indicating that any pleiotropic bias is likely minor.

### Heterogeneity Statistics


::: {.cell}

```{.r .cell-code}
# Heterogeneity tests for IVW and MR-Egger
heterogeneity <- mr_heterogeneity(data_steiger)
heterogeneity|>
  select(-starts_with('id')) |> 
  mutate(
    I2 = pmax(0, (Q - Q_df) / Q) * 100 # 
  ) |>
  kable(caption="MR Heterogeneity Results for Calcium Effects on Cholesterol",
        digits=c(0,0,0,3,3,99))
```

::: {.cell-output-display}


Table: MR Heterogeneity Results for Calcium Effects on Cholesterol

|outcome                            |exposure             |method                    |       Q| Q_df|       Q_pval| I2|
|:----------------------------------|:--------------------|:-------------------------|-------:|----:|------------:|--:|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |MR Egger                  | 604.623|  273| 3.710043e-27| 55|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |Inverse variance weighted | 606.474|  274| 3.320526e-27| 55|


:::

```{.r .cell-code}
# Columns: method, Q, Q_df, Q_pval
# Interpretation: small Q_pval (<0.05) indicates heterogeneity among SNPs
```
:::


This is expected with polygenic traits and does not necessarily invalidate the overall causal estimate, particularly since robust methods (weighted median, weighted mode) gave consistent results.


::: {.cell}

```{.r .cell-code}
single_snp_results <- mr_singlesnp(data_steiger)
mr_funnel_plot(single_snp_results) 
```

::: {.cell-output .cell-output-stdout}

```
$`Calcium (UK Biobank).LDL Cholesterol (MGI-BioVU LabWAS)`
```


:::

::: {.cell-output-display}
![](figures/calcium-ldlc-funnel-plot-1.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```

attr(,"split_type")
[1] "data.frame"
attr(,"split_labels")
           id.exposure                         id.outcome
1 Calcium (UK Biobank) LDL Cholesterol (MGI-BioVU LabWAS)
```


:::

```{.r .cell-code}
# Get overall IVW estimate for the vertical line
ivw_beta <- calcium.ldlc.mr |> filter(method=="Inverse variance weighted (multiplicative random effects)") |> pull(b)

# Determine y-range based on your data
y_min <- 0
y_max <- max(single_snp_results$se^{-1}) * 1.1  # 10% padding above max precision

# Generate a fine grid of precision values
precision_grid <- seq(y_min, y_max, length.out = 1000)

# Compute boundaries: ivw_beta ± 1.96 / precision
lower_bound <- ivw_beta - 1.96 / precision_grid
upper_bound <- ivw_beta + 1.96 / precision_grid

# Create data frame for boundaries
bounds_df <- data.frame(precision = precision_grid, lower = lower_bound, upper = upper_bound)

# Plot
ggplot(single_snp_results, aes(x = b, y = 1/se)) +
  # Scatter points for each SNP
  geom_point(size = 1) +
  # Vertical line at IVW estimate
  geom_vline(xintercept = ivw_beta, linetype = "solid", color = "#ff7f0e", size = 1) +
  # Curved pseudo-95% CI boundaries (the cone)
  geom_line(data = bounds_df, aes(x = lower, y = precision), linetype = "dashed") +
  geom_line(data = bounds_df, aes(x = upper, y = precision), linetype = "dashed") +
  # Customize axes and labels
  labs(
    x = "Estimate (Beta-IVW)",
    y = "Precision (1/Standard Error)",
    title = ""
  ) +
  # Apply clean theme and limit y to >=0
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, y_max), xlim = c(min(single_snp_results$b), max(single_snp_results$b)))  # Adjust x-limits for visibility
```

::: {.cell-output-display}
![](figures/calcium-ldlc-funnel-plot-2.png){width=672}
:::
:::


### Leave-one-out Analysis

Using IVW methods


::: {.cell}

```{.r .cell-code}
# LOO using IVW
loo_res <- mr_leaveoneout(data_steiger)
loo_res |> 
  mutate(diff = b - filter(calcium.ldlc.mr, method=="Inverse variance weighted (multiplicative random effects)")$b) |>
  arrange(-abs(diff)) |>
  head() |>
  select(SNP,diff,b,se,p) |>
  kable(caption="Leave-One-Out Results for Calcium-LDL Cholesterol Analysis (IVW method) for Influential SNPs",
        digits=c(0,5,5,5,5))
```

::: {.cell-output-display}


Table: Leave-One-Out Results for Calcium-LDL Cholesterol Analysis (IVW method) for Influential SNPs

|SNP              |     diff|        b|      se|       p|
|:----------------|--------:|--------:|-------:|-------:|
|1:109817192:A:G  | -0.01386| -0.01495| 0.03198| 0.64027|
|2:27742603:T:C   | -0.01153| -0.01261| 0.03847| 0.74298|
|8:9183358:A:G    | -0.00957| -0.01065| 0.03834| 0.78112|
|11:126218773:C:G |  0.00940|  0.00832| 0.03811| 0.82722|
|3:121993247:A:G  |  0.00907|  0.00798| 0.04083| 0.84498|
|19:19717056:A:G  |  0.00554|  0.00446| 0.03792| 0.90646|


:::

```{.r .cell-code}
# Columns: SNP, nsnp, b, se, pval — gives causal estimate with each SNP removed once

# Optional: plot LOO results

ggplot(loo_res, aes(x = reorder(SNP, -b), y = b)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin = b - 1.96*se, ymax = b + 1.96*se), width = 0.01 ,alpha=0.5) +
  coord_flip() +
  labs(x = "SNP Removed", y = "Estimate (Beta-IVW; leave-one-out)") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic(base_size=16) +
  theme(axis.text.y = element_text(size = 1)) 
```

::: {.cell-output-display}
![](figures/calcium-ldlc-mr-loo-1.png){width=672}
:::
:::

Leave-one-out analyses suggested that four SNPs had a relatively large influence on the IVW estimate, but removal of either SNP did not qualitatively change the overall conclusion.


### MR-CAUSE Analysis

CAUSE was used to model both correlated and uncorrelated horizontal pleiotropy.  Correlated pleiotropy are the effects of the SNPs an outcome not through the trait but through a confounder.  Uncorrelated horizontal pleiotropy is direct effects of the SNPs on the outcome independent of the modeled trait.  This is described in [@morrisonMendelianRandomizationAccounting2020].


::: {.cell}

```{.r .cell-code}
#devtools::install_github("jean997/cause@v1.2.0")
library(cause)
tc.cause.data <-
  data_steiger |>
  rename(
    snp = SNP,
    beta_hat_1 = beta.exposure,
    beta_hat_2 = beta.outcome,
    seb1 = se.exposure,
    seb2 = se.outcome
  ) |>
  new_cause_data()

tc.params_ests <- est_cause_params(
  X = tc.cause.data,                    # Merged data
  variants = tc.cause.data$snp,
  optmethod = "mixSQP",     # Default & recommended
  null_wt = 10,             # Weight on null (default)
  max_candidates = Inf      # Full grid (default)
)
```

::: {.cell-output .cell-output-stdout}

```
Estimating CAUSE parameters with  277  variants.
1 0.1077897 
2 0.001350809 
3 1.888131e-05 
4 1.241148e-07 
5 1.009747e-07 
6 8.446399e-08 
```


:::

```{.r .cell-code}
calcium.ldlc.cause <- cause(X=tc.cause.data,
                               param_ests = tc.params_ests)
```

::: {.cell-output .cell-output-stdout}

```
Estimating CAUSE posteriors using  277  variants.
```


:::

```{.r .cell-code}
calcium.ldlc.cause$elpd |> kable(caption="if delta_elpd is negative, model2 is a better fit, in this case means the causal model is better than the pleiotropic sharing model or either null models")
```

::: {.cell-output-display}


Table: if delta_elpd is negative, model2 is a better fit, in this case means the causal model is better than the pleiotropic sharing model or either null models

|model1  |model2  | delta_elpd| se_delta_elpd|         z|
|:-------|:-------|----------:|-------------:|---------:|
|null    |sharing |  0.3723442|     0.6675411| 0.5577846|
|null    |causal  |  1.2815556|     0.9909655| 1.2932394|
|sharing |causal  |  0.9092115|     0.6428994| 1.4142361|


:::

```{.r .cell-code}
plot(calcium.ldlc.cause, type="data",intern=TRUE) -> tmp.plots
tmp.plots[[1]]
```

::: {.cell-output-display}
![](figures/calcium-ldlc-cause-1.png){width=672}
:::

```{.r .cell-code}
tmp.plots[[2]]
```

::: {.cell-output-display}
![](figures/calcium-ldlc-cause-2.png){width=672}
:::

```{.r .cell-code}
tmp.plots[[3]]
```

::: {.cell-output-display}
![](figures/calcium-ldlc-cause-3.png){width=672}
:::

```{.r .cell-code}
summary(calcium.ldlc.cause, ci_size = 0.95)$tab |> kable(caption="Pathway estimates and 95% confidence interveals for estimated effect sizes, ")
```

::: {.cell-output-display}


Table: Pathway estimates and 95% confidence interveals for estimated effect sizes, 

|model   |gamma              |eta                 |q              |
|:-------|:------------------|:-------------------|:--------------|
|Sharing |NA                 |-0.58 (-2.27, 0.8)  |0.03 (0, 0.18) |
|Causal  |0.01 (-0.08, 0.09) |-0.57 (-1.76, 0.75) |0.05 (0, 0.24) |


:::
:::


From the CAUSE analyses there is qualitative evidence to prefer the causal pathway compared with the shared (pleiotropic) pathways (p=0.9213537). The estimated causal effect ($\gamma$) is 0.01 (-0.08, 0.09) and the residual correlated pleiotropy was minimal after accounting for this causal effect. The $\eta$ = -0.57 (-1.76, 0.75) is near zero for the causal model but is slightly larger for the sharing model [$\eta$=-0.58 (-2.27, 0.8)]. To explain this data without a causal effect, CAUSE would require more correlated pleiotropy.  In the absence of a causal effect (sharing model), correlated horizontal pleiotropy would explain 0.05 (0, 0.24)% of the SNPs would require correlated pleiotropy for the causal model, but 0.03 (0, 0.18)% of the SNPs would. 

Alternate explanation with assistance from ChatGPT:

The CAUSE model comparison favored the causal model over both the null and sharing models, although none of the differences reached statistical significance. For example, comparing the sharing vs. causal models yielded a $\Delta$ELPD (Expected Log Pointwise Predictive Density) of 0.9092115 with a standard error of 0.6428994 (z score of = 1.4142361). The causal model estimated a null effect of calcium on LDL-C (0.01 (-0.08, 0.09)), while the corresponding pleiotropic parameter $\eta$ was centered near zero (-0.57 (-1.76, 0.75)), suggesting minimal directional pleiotropy. The estimated fraction of variants exhibiting correlated pleiotropy (q) was small under the causal model (0.05 (0, 0.24)), and lower than under the sharing model (0.03 (0, 0.18)). 

### Summary of Analyses


::: {.cell}

```{.r .cell-code}
calcium.ldlc.cause.summary <- 
  summary(calcium.ldlc.cause, ci_size = 0.95)$tab |> 
  as_tibble() |>
  filter(model=="Causal") |>
  mutate(method=fct_recode(as.factor(model), "MR-CAUSE"="Causal")) |>
  select(method,gamma) |>
  separate(
    col = gamma, into = c("b", "ci"), sep = " \\(",remove = TRUE) |>
  mutate(
    ci = str_remove(ci, "\\)$"),          # remove trailing ")"
    ci = str_squish(ci)) |>              # clean any extra spaces
  separate(ci, into=c("lower.ci","upper.ci"), sep=", ") |>
  mutate(se = (as.numeric(upper.ci)-as.numeric(lower.ci))/2/1.96) |>
  mutate(b=summary(calcium.ldlc.cause)$quants[[2]][1,'gamma']) |>
  select(method,b,se)

method.order <- c("IVW-RE",
                  "IVW-FE",
                  "Weighted median",
                  "MR Egger",
                  "Weighted mode",
                  "MR-PRESSO (Raw)",
                  "MR-PRESSO (Corrected)",
                  "MR-RAPS",
                  "MR-CAUSE")

calcium.ldlc.summary <-
  calcium.ldlc.mr.mrpresso |> 
  select(-starts_with('id')) |>
  bind_rows(calcium.ldlc.cause.summary) |>
  mutate(method=fct_recode(as.factor(method),
                                  "IVW-RE"="Inverse variance weighted (multiplicative random effects)",
                                  "IVW-FE"="Inverse variance weighted (fixed effects)",
                           "MR-RAPS"="Robust adjusted profile score (RAPS)",
                           "MR-PRESSO (Corrected)" = "MR-PRESSO (Outlier-corrected)")) |>
  mutate(method = factor(method, levels=method.order)) |>
  arrange(method) |>
  fill(outcome,exposure,nsnp) 
  
  
calcium.ldlc.summary |>   
  kable(caption="MR Results for Calcium on LDL-C",
        digits=c(0,0,0,0,4,4,99))
```

::: {.cell-output-display}


Table: MR Results for Calcium on LDL-C

|outcome                            |exposure             |method                | nsnp|       b|     se|      pval|
|:----------------------------------|:--------------------|:---------------------|----:|-------:|------:|---------:|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |IVW-RE                |  275| -0.0011| 0.0383| 0.9773529|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |IVW-FE                |  275| -0.0011| 0.0258| 0.9663122|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |Weighted median       |  275| -0.0599| 0.0484| 0.2162812|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |MR Egger              |  275| -0.0659| 0.0806| 0.4141755|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |Weighted mode         |  275| -0.0004| 0.0636| 0.9951880|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |MR-PRESSO (Raw)       |  275|  0.0007| 0.0383| 0.9858026|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |MR-PRESSO (Corrected) |  275| -0.0016| 0.0295| 0.9557368|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |MR-RAPS               |  275| -0.0044| 0.0328| 0.8927706|
|LDL Cholesterol (MGI-BioVU LabWAS) |Calcium (UK Biobank) |MR-CAUSE              |  275|  0.0115| 0.0434|        NA|


:::

```{.r .cell-code}
calcium.ldlc.summary |> 
  write_csv("MR Results - Calcium - LDL Cholesterol.csv")

calcium.ldlc.summary |>
  mutate(method = factor(method, levels = rev(method.order))) %>% #reverse order
  ggplot(aes(y=method ,x=b)) +
  geom_point() +
  geom_errorbar(aes(xmin=b-1.96*se, xmax=b+1.96*se), width=0.2) +
  theme_classic(base_size=16) +
  labs(title="",
       y="",
       x="Effect Size (Beta)") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") 
```

::: {.cell-output-display}
![](figures/ldlc-calcium-mr-summary-1.png){width=672}
:::
:::



## Session Information


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
 [1] cause_1.2.0        MRPRESSO_1.0       ggrepel_0.9.6      TwoSampleMR_0.6.22
 [5] knitr_1.50         lubridate_1.9.4    forcats_1.0.1      stringr_1.6.0     
 [9] dplyr_1.1.4        purrr_1.2.0        readr_2.1.6        tidyr_1.3.1       
[13] tibble_3.3.0       ggplot2_4.0.1      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] gtable_0.3.6          xfun_0.54             htmlwidgets_1.6.4    
 [4] psych_2.5.6           lattice_0.22-7        tzdb_0.5.0           
 [7] vctrs_0.6.5           tools_4.5.2           generics_0.1.4       
[10] curl_7.0.0            parallel_4.5.2        pkgconfig_2.0.3      
[13] Matrix_1.7-4          SQUAREM_2021.1        data.table_1.17.8    
[16] RColorBrewer_1.1-3    S7_0.2.1              RcppParallel_5.1.11-1
[19] truncnorm_1.0-9       lifecycle_1.0.4       rootSolve_1.8.2.4    
[22] compiler_4.5.2        farver_2.1.2          mnormt_2.1.1         
[25] htmltools_0.5.8.1     mr.raps_0.4.2         yaml_2.3.10          
[28] pillar_1.11.1         crayon_1.5.3          nlme_3.1-168         
[31] rsnps_0.6.1           tidyselect_1.2.1      digest_0.6.38        
[34] nortest_1.0-4         stringi_1.8.7         ashr_2.2-63          
[37] labeling_0.4.3        splines_4.5.2         fastmap_1.2.0        
[40] grid_4.5.2            invgamma_1.2          cli_3.6.5            
[43] magrittr_2.0.4        loo_2.8.0             crul_1.6.0           
[46] withr_3.0.2           scales_1.4.0          bit64_4.6.0-1        
[49] timechange_0.3.0      rmarkdown_2.30        matrixStats_1.5.0    
[52] bit_4.6.0             gridExtra_2.3         hms_1.1.4            
[55] evaluate_1.0.5        irlba_2.3.5.1         mgcv_1.9-4           
[58] rlang_1.1.6           mixsqp_0.3-54         Rcpp_1.1.0           
[61] glue_1.8.0            httpcode_0.3.0        rstudioapi_0.17.1    
[64] vroom_1.6.6           jsonlite_2.0.0        R6_2.6.1             
[67] plyr_1.8.9            intervals_0.15.5     
```


:::
:::

