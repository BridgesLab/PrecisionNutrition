## Data Direct Data Pull

## Filters

- Initial: 5,529,222 patients
- Michigan Medicine Primary Care and MGI: 111,011 (this is mostly to get the download to <100,000 participants)
- 18+ at Encounter (2000-2026): 109,820
- Has some cholesterol or triglyceride measurement: 51,122
- Also has some calcium measurement: 50,225

### Data

- Demographics
   - Gender
   - Race
   - Ethnicity
   - DoB 
   - Death (can be validated to death registry)
   - Cause of Death (can be validated to death registry)
- Social History
   - Drinking status
   - Smoking status
- Lab Results
   - Albumin, Calcium, Cholesterol, Triglycerides, HDL, LDL (variety of names), Vitamin D, Parathyroid Hormone, B-CTX, NTX, PINP, osteocalcin, phosphate, phosphorous, creatinine, eGFR
   - Value
   - Units
   - Range
   - Collection Date
- Comorbidities
   - Charlson - Congestive Heart Failure, Diabetes (2x), Renal Disease, Myocardial Infarction, Peripheral Vascular Disease, Total Score, Age Adjusted Total Score
   - Elixhauser - Alcohol Abuse, CardiacArrythmia, Congestive Heart Failure, Diabetes (2x), Hypertension (2x), Obesity, Pulmonary Circulation Disorders, Renal Failure, Total Score
- EncounterAll - Admit Date, AgeInYears
- Encounter Anthropometrics
   - BMI
   - Height
   - Weight
- Medication
   - Statins
   - Date, MinDose, MaxDose, UMDose, Frequency, Status, Strength
- GISNeighborhood Affluence
   - affluence
   - disadvantage
   - education
- NursingVitalSigns
   - MAP, SBP, DBP, ObservationDate
- Diagnosis Comprehensive
   - MACE (see [MACE_ICD_Codes_List.md](MACE_ICD_Codes_List))
      - ICD9: 410 + (433 + 434 + 436) + (430 + 431 + 432) 
      - ICD-10: I21 + I22 + I45 + I63 + I64 + I60 + I61 + I62 + I46 (no I64)
   - Osteoporosis and Related Fractures [1]
      - ICD‑10: M80.x, M81.x
      - ICD‑9 733.0x (733.00 unspecified; 733.01 senile; 733.02 idiopathic; 733.03 disuse; 733.09 other).

[1] Chen Y, Harrold LR, Yood RA, Field TS, Briesacher BA. Identifying patients with osteoporosis or at risk for osteoporotic fractures. Am J Manag Care. 2012 Feb 1;18(2):e61-7. PMID: 22435886; PMCID: PMC4841251.