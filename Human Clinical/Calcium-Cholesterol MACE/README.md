## Analysis of Data for Cholesterol and MACE/Osteoporosis Outcomes

- Initial: 5,529,222 patients
- Michigan Medicine at Primary Care: 352,177
- Has at least one LDL measurement since 1/1/2000 201,073

This was going to need to be done in multiple batches by parsing patient ID's and then resubmitting this query in multiple rounds (excluding the prior ones).  First query was males (91,653), second was females/unknown - married (50,042) and third was females/unknown - not married/unknown (59,378 patients).  

### Data Cleaning

1.  The first step was to combine the three files.  This was done in the `combining_raw_data.qmd` script.  This generated new combined files in the `combined_data` folder with the same filenames.
2.  Diagnoses are cleaned via the `diagnosis_cleaning.qmd` script.  It generates a new file in `combined_data` called **DiagnosesCleaned.csv**.  It uses the DiagnosesComprehensiveAll.csv, EncounterAll.csv, and EncounterAnthropometricsBMI.csv files to annotate this file.  It generates the MACE and Osteoporosis diagnoses.
3.  Lab results are cleaned via the `labs_cleaning.qmd` script.  This generates a new file in `combined_data` called **LabResultsCleaned.csv**.  It uses the LabResults.csv, EncounterAll.csv, and EncounterAnthropometricsBMI.csv files to annotate this file.  It generates the MACE and Osteoporosis diagnoses.  This script combines similar measurements into an aggregated column called `test_name`, and adjusts units to mg/dL when not already.
4.  Statin data are cleaned via the `medication_cleaning.qmd` script.  This generates a new file in `r combined_data` called **MedicationOrdersCleanedStatins.csv**.  It only uses the MedicationOrdersComprehensive.csv file.  It generates statin intervals for use as a time-varying covariate in Cox proportional hazard models.

The data had to be pulled int