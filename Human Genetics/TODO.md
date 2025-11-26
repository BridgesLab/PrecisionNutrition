## TODO For Manuscript

### Missing Data

- [ ] n for MGI-BioVU (add to scripts/manuscript/Figure 1A)
- [x] population details
- [ ] Units for beta coefficients from MGI-BioVU (SD for UKBB)

### Mising Analsyes

- [ ] Regenerate SNP lists for Total Cholesterol on Fracture Risk and BMD
- [ ] Add MR-CAUSE analyses for models to test for correlated pleiotropy
- [ ] Add MR-PRESSO analyses 
- [ ] Add MR-RAPS analyses 
- [ ] Use phenoscanner to identify potential confounders


## Limtations Notes

### MR Assumptions

- **Relevance assumption**: that genetic variants are strongly associated with the exposure of interest, test 
- **Independence assumption**: variants should not be associated with confounders of the relatoinship.  Test using `phenoscanner` or similar tools
- **Exclusion restriction assumption**: Test that MR-Egger intercept tests do not detect directional pleiotropy.

### Pleiotropy 

SNPs affecting two or more seemingly unrelated phenotypes, because genes can participate in multiple biological pathways

#### Vertical Plieotropy

Variant is related to the exposure, this is what we are evaluating.

#### Directional Uncorrelated Horizontal Pleiotropy

SNPs have effects independent of the effects on the exposure.  Creates an average directional effect on the outcome.  Assess by Egger's intercept != 0 (if the InSIDE assumption holds) as well as MR-PRESSO Global Test.  As a general rule, if the global test is significant and the distortion test is not, we should use the MR-PRESSO corrected estimate.  If the global test if not significant we should use the IVW estimate.  

##### InSIDE Assumption

The InSIDE assumption (that instrument strength is independent of any direct (pleiotropic) effect on the outcome) is a factor for MR-Egger but not MR-PRESSO or IVW methods.  As a general rule the stronger a SNP is as an instrument for the exposure, the less (or equally) likely it is to have a direct effect on the outcome that bypasses the exposure.  If using MR-Egger this can be tested via the Sanderson-Windmeijer Conditional F-statistic.  A value <10 indicates it is probably violated.

#### Correlated (Coordinated) Pleiotropy

SNPS influence the exposure and outcome through a shared mechanism or confounder.  Can produce a false positive causal effect.  Can fit using CAUSE models (Causal Analysis Using Summary Effect Estimates https://github.com/jean997/cause)

#### Balanced Horizontal Pleiotropy

Genes have multidirectional effects, but they cancel each other out.  Increases heterogeneity, but no directional bias.  Can use MR-PRESSO (Pleiotropy RESidual Sum and Outlier) analyses to detect outlier driven and balanced pleiotropy.  A significant global test, or multiple outliers supports average pleiotropy.  Also if MR-RAPS causal estiamte remain stable relative to IVW, pleiotropy is likely balanced rather than directional.

#### Assumptions of Models

| Method | Assumptions |
|--------|----------------------|
| IVW |  Assumes all variants are valid IVs or that any pleiotropic effects are balanced (mean zero); directional pleiotropy biases IVW. |
| MR-Egger | Allows all variants to be invalid, but assumes uncorrelated (InSIDE) pleiotropy and tests for directional average pleiotropy via the intercept |
| Weighted Median | Gives a consistent estimate if ≥50% of the weight comes from valid instruments and no single SNP dominates the weight. |
| Weighted Mode | Consistent if the largest “cluster” (plurality) of SNPs, by weight, shares the same causal effect and is valid. |
| MR-PRESSO | Detects global horizontal pleiotropy, identifies outlier SNPs driving it, and provides an outlier‑corrected IVW estimate plus a distortion test. |
| MR-RAPS | Designed for many weak instruments and robust to some pleiotropy via a robust loss; also useful when there is substantial measurement error in SNP–exposure effects. |
| MR-CAUSE | Explicitly models both uncorrelated and correlated pleiotropy and compares a causal vs sharing (pleiotropy) model using genome‑wide summary data |

#### Decision Tree for Model Selection

| Step | Process | Decision / Threshold | Action |
|---|---|---|---|
| **1. Instrument Validation** | **Analyze Genetic Instruments** | **Mean F-statistic ≥ 10** | **Strong Instruments.** Proceed to MR. |
| | | **Mean F-statistic < 10** | **Weak Instruments.** Use **MR-RAPS** (robust to weak instruments) as primary. Interpret all results cautiously. |
| **2. Primary Analysis & Heterogeneity** | **Fit IVW Model & Assess Heterogeneity** (Cochran's Q, I²) | **Q p ≥ 0.05 AND I² < 50-60%** | **No significant heterogeneity.** Use **IVW (Fixed Effects, FE)** as the primary estimate. |
| | | **Q p < 0.05 OR I² > 50-60%** | **Heterogeneity present.** Suspect pleiotropy/outliers. Use **IVW (Random Effects, RE)** as the primary estimate. |
| **3. Directional Pleiotropy** | **Check MR-Egger Intercept** | **Intercept p ≥ 0.05** (CI includes 0) | **No significant directional pleiotropy.** IVW-RE/FE is valid. Report Weighted Median/Mode and Egger slope as sensitivity analyses. |
| | | **Intercept p < 0.05** (CI excludes 0) | **Directional pleiotropy likely.** Emphasize **Weighted Median** and **MR-Egger slope** as the most robust estimates. |
| **4. Outlier-Driven Pleiotropy** | **Run MR-PRESSO** | **Global Test p ≥ 0.05** | **No significant outliers.** Prior choice of IVW/Median/Mode stands. |
| | | **Global Test p < 0.05** (Outliers Present) | **Outliers Detected.** Proceed to Distortion Test. |
| | | **Distortion Test p < 0.05** (Estimate Changed) | Use **MR-PRESSO Outlier-Corrected IVW** as the primary estimate, alongside Median/Mode. |
| | | **Distortion Test p ≥ 0.05** (Estimate Unchanged) | **IVW-RE** is still acceptable, but acknowledge and report sensitivity due to outliers. |
| **5. Sensitivity to Correlated Pleiotropy** | **Check for Correlated Pleiotropy** (e.g., genetic correlation, CAUSE) | **CAUSE "Causal" Model** (or no strong correlation) | **Final interpretation:** Interpret the most robust/primary MR estimates selected from Steps 2-4. |
| | | **CAUSE "Sharing" Model** favored | **Reframed Interpretation:** Evidence for **shared genetic architecture**, not clear causality. Interpret results cautiously. |

## Other Notes
- IVW random-effects will serve as the **primary estimator**.
- All other MR methods reported for sensitivity.
- Same harmonisation and filtering pipeline applied in both directions.
- Pleiotropy and heterogeneity results should be cross-referenced in main text as quality checks.
- Positive control (Calcium→Calcium) included in Supplement for validation.

## Potential Downstream Analyses
- Fracture Risk -- https://pheweb.org/MGI/pheno/X800.1, https://pheweb.org/MGI/pheno/X743.22, https://pheweb.org/MGI/pheno/X743.21
- BMD (total, lumbar spine, femoral neck) -- not in MGI-BioVU or MGI
- Osteoporosis (via dexa or ICD) -- https://pheweb.org/MGI/pheno/X743
- Urinary calcium excretion/Nephrolithiasis -- not in MGI-BioVU or MGI
- CKD/eGFR -- https://pheweb.org/MGI-BioVU/pheno/Creat
- Bone turnover markers (osteocalcin, BSAP, P1NP CTX or NTX telopeptides) -- not in MGI-BioVU
- Cognitive decline/Alzheimer's disease -- https://pheweb.org/MGI/pheno/X290.11 (only a weak rationale)
- Vacular calficication, coronary artery calcium score, heart disease -- too much horizontal pleiotropy
- Iscemic > hemmorrhagic stroke


### Mechanistic insights
 - 25-hydroxyvitamin D -- https://pheweb.org/MGI-BioVU/pheno/Vit-D
 - FGF23 -- not in MGI-BioVU
 - PTH -- not in MGI-BioVU
 - Estrogen, SHBG, testosterone -- not in MGI-BioVU
