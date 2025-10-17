#

## Missing Data

- [ ] n for MGI-BioVU (add to scripts/manuscript/Figure 1A)
- [x] population details
- [ ] Units for beta coefficients from MGI-BioVU (SD for UKBB)

## Analysis Pipeline

- [ ] Instrument QC table (pre)
- [ ] Expand SP2 & Flow chart counts
- [ ] Harmonise (produce harmonised table) → write harmonised SNP table (step 4)
- [ ] Recompute per-SNP R² & F on harmonised set and write post-harmonisation QC table (step 4)
- [ ] Scatter plot (step 5)
- [ ] MR estimates table + forest (steps 6–7)
- [ ] Heterogeneity + pleiotropy + funnel (steps 8–9)
- [ ] LOO + radial / MR-PRESSO if needed (steps 10–11)
- [ ] Steiger summary + proxies table (steps 12–13)
- [ ] Bidirectional summary figure/table (step 14)
- [ ] Write supplement SNP table (step 15)

## Proposed Manuscript Structure

Generated via ChatGPT

This plan summarizes the figures and tables to include in the main manuscript and supplementary materials for the bidirectional MR analysis of calcium and cholesterol (LDL and total cholesterol).  
It follows STROBE-MR reporting principles and emphasizes clarity, reproducibility, and visual balance between main and supplementary materials.

### **Figure 1. Study Design and Instrument Quality**
**Panel A – Study schematic**
- Diagram showing bidirectional MR design:
  - Calcium → LDL & TC
  - LDL & TC → Calcium
- Include positive control analysis (Calcium UK Biobank → Calcium Michigan)
- Indicate data sources, sample sizes, and harmonization steps.

**Panel B – Instrument summary**
- Compact table (or inset) summarizing instrument strength:
  | Exposure | N SNPs | Cumulative R² | Mean F | Median F | Mean MAF |
  |-----------|--------|---------------|---------|-----------|-----------|
  | Calcium (UKBB) |   |   |   |   |   |
  | LDL |   |   |   |   |   |
  | Total Cholesterol |   |   |   |   |   |

> **Purpose:** Orient reader; demonstrate adequate instrument strength and data independence.

---

### **Figure 2. Primary Bidirectional MR Results**
Four panels arranged as a 2×2 grid:

|               | LDL Outcome | TC Outcome |
|----------------|--------------|-------------|
| **Calcium Exposure** | (A) Forest plot: IVW, Egger, Weighted Median, Mode | (B) Forest plot: same methods |
| **Reverse Direction (LDL/TC Exposure)** | (C) Forest plot: LDL → Calcium | (D) Forest plot: TC → Calcium |

Each panel:
- Shows β estimates ± 95% CI for each MR method.
- Highlights IVW as the primary estimate (bold or color-coded).

> **Purpose:** Compare causal direction and consistency across MR estimators.

---

### **Figure 3. Sensitivity and Robustness Diagnostics**
Three panels for the primary direction (Calcium → LDL):

**Panel A – Scatter plot**
- SNP-level β_exposure vs β_outcome.
- Overlay IVW and MR-Egger regression lines.

**Panel B – Leave-One-Out (LOO) influence plot**
- Each SNP omitted in turn; show shifts in IVW estimate.

**Panel C – Funnel plot**
- Symmetry check for heterogeneity / directional pleiotropy.

> **Purpose:** Demonstrate overall robustness, absence of single-SNP dominance, and evaluate pleiotropy visually.

---

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
