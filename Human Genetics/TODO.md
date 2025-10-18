#

## Missing Data

- [ ] n for MGI-BioVU (add to scripts/manuscript/Figure 1A)
- [x] population details
- [ ] Units for beta coefficients from MGI-BioVU (SD for UKBB)

## Proposed Manuscript Structure

Generated via ChatGPT

This plan summarizes the figures and tables to include in the main manuscript and supplementary materials for the bidirectional MR analysis of calcium and cholesterol (LDL and total cholesterol).  
It follows STROBE-MR reporting principles and emphasizes clarity, reproducibility, and visual balance between main and supplementary materials.

### **Table 1. Summary of Main MR Results**
Compact table of key causal estimates.

| Direction | Exposure | Outcome | nsnp | Method | Beta | SE | 95% CI | p | Steiger direction |
|------------|-----------|----------|------|---------|------|----|--------|---|-------------------|
| 1 | Calcium | LDL | | IVW | | | | | |
| 2 | Calcium | TC | | IVW | | | | | |
| 3 | LDL | Calcium | | IVW | | | | | |
| 4 | TC | Calcium | | IVW | | | | | |

> **Purpose:** Single summary table showing directionality, effect sizes, and significance.

---

## SUPPLEMENTARY MATERIALS

### **Supplementary Figure S1. Instrument selection flow diagram**
- Summary of SNP selection and filtering

---

### **Supplementary Figure S2. Positive Control Analysis**
Confirms validity of analysis pipeline:
-  A. Scatter plot showing beta coefficients
-  B. Forest plots showing Calcium (UKBB) → Calcium (Michigan).
-  C. Leave-one-out plot: usually unnecessary for a positive control unless a reviewer asks.
-  D. Funnel plot / heterogeneity check: can be skipped unless you want to emphasize no pleiotropy.

---

### **Supplementary Figure S3. Full Scatter and Funnel Plots**
- Scatter and funnel plots for **all four** analyses:
  1. Calcium → LDL  
  2. Calcium → TC  
  3. LDL → Calcium  
  4. TC → Calcium

---

### **Supplementary Figure S4. MR-PRESSO and Radial MR Analyses**
- Radial plots identifying outlier SNPs and post-correction results.
- Optional if using MR-PRESSO or radialMR package.

---

### **Supplementary Figure S5. Heterogeneity and Pleiotropy Checks**
- Barplots or small tables of:
  - Cochran’s Q (IVW and Egger)
  - MR-Egger intercept, SE, p-value

---

### **Supplementary Figure S6. Influential SNPs and Annotations**
- Leave-one-out tables/plots identifying top 10 influential SNPs per analysis.
- Include nearest gene annotation for context.

---

### **Supplementary Table S1. Full Instrument List**
- SNP-level data for each exposure:
  - CHR, POS, REF, ALT, effect_allele, other_allele,
    beta.exposure, se.exposure, p.exposure, eaf.exposure, R², F, clump_lead.

---

### **Supplementary Table S2. Harmonisation Summary**
- Counts of SNPs at each filtering stage:
  - Candidate → Harmonised → Palindromic dropped → Steiger filtered → Final N.

---

### **Supplementary Table S3. Complete MR Output**
- Full results for all MR methods:
  - IVW, MR-Egger, Weighted Median, Weighted Mode, MR-RAPS, MR-PRESSO (if used)
  - For each direction and outcome.

---

### **Supplementary Methods**
- Detailed pipeline description:
  - Clumping parameters (r², window)
  - Reference panel
  - Harmonisation strategy
  - Proxy SNP handling
  - Software and package versions
  - Exact GWAS summary statistic sources and sample sizes
  - Code availability (GitHub/Zenodo)

---

## Notes
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
