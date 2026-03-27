## Data Direct Data Pull

## Filters

### Calcium and Cholesterol

- Initial: 5,529,222 patients
- MGI: 111,011 (this is mostly to get the download to <100,000 participants)
- 18+ at Encounter (2000-2026): 109,820
- Has some cholesterol or triglyceride measurement: 51,122
- Also has some calcium measurement: 50,225

### Cholesterol and Outcomes

- Initial: 5,529,222 patients
- Michigan Medicine at Primary Care: 352,177
- Has at least one LDL measurement since 1/1/2000 201,073

This was going to need to be done in multiple batches by parsing patient ID's and then resubmitting this query in multiple rounds (excluding the prior ones).  First query was males (91,653), second was females/unknown - married (50,042) and third was females/unknown - not married/unknown (59,378 patients).  These were combined prior to all scripts.

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