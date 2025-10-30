## TODO For Manuscript

### Missing Data

- [ ] n for MGI-BioVU (add to scripts/manuscript/Figure 1A)
- [x] population details
- [ ] Units for beta coefficients from MGI-BioVU (SD for UKBB)

### Mising Analsyes

- [ ] Regenerate SNP lists for Total Cholesterol on Fracture Risk and BMD
- [ ] Add CAUSE analyses for models to test for correlated pleiotropy
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

| Method | Pleiotropy Assumption |
|--------|----------------------|
| IVW | No pleiotropy |
| MR-Egger | Uncorrelated pleiotropy |
| Weighted Median | <50% invalid instruments |
| Weighted Mode | Plurality of instruments valid |
| MR-PRESSO | Outlier correction |
| MR-RAPS | Weak instruments |
| CAUSE | Correlated pleiotropy |

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
