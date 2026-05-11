# Analysis of Summary Statistics for Calcium and LDL-C from UK Biobank

We wanted to be able to download summary statistics for some traits, and do clumping analysis using plink2.

## Reproducibility Notes

The LD reference panel and clumped instrument files are stored on GreatLakes/DataDen and are not part of this repository. The MGI-BioVU LabWAS summary statistics (`PheWeb Summary Statistics/`) are restricted to University of Michigan researchers and cannot be shared publicly; external users would need to substitute comparable publicly available GWAS summary statistics. Exact software versions for bcftools, samtools, and tabix are not pinned in the SLURM scripts (loaded via the GreatLakes module system); check `sessionInfo()` output in each rendered `.qmd` for the R environment used at time of analysis.

## Script Execution Order

The analyses must be run in the following order:

1. `reference-build.sh` (via `reference-build.slurm`) — builds 1000G EUR LD reference panel on HPC
2. `ld-clumping.slurm` — downloads Pan-UKBB summary statistics and performs LD clumping against the reference panel
3. `ld_clump_analysis.qmd` — filters instruments by MAF and writes out final SNP lists
4. Main MR analyses (independent, can be run in any order):
   - `mr-tc-calcium.qmd`
   - `mr-ldlc-calcium.qmd`
   - `mr-calcium-tc.qmd`
   - `mr-calcium-ldlc.qmd`
5. `drug_target_mr_analysis.qmd` — requires internet access to OpenGWAS API
6. `drug_target_mr_ukb.qmd` — requires internet access to OpenGWAS API
7. `mr-downstream-analyses.qmd` — requires internet access to OpenGWAS API
8. `summary-tables.qmd` — collates outputs from all MR scripts above

## Software Requirements

* plink2 (v2.00a5.12LM AVX2 AMD (25 Jun 2024))
* bcftools (version not pinned; loaded via GreatLakes module system)
* samtools (version not pinned; loaded via GreatLakes module system)
* tabix (version not pinned; loaded via GreatLakes module system)

### R Packages

The following R packages are required for the `.qmd` analysis scripts.  CRAN packages can be installed with `install.packages()`; Bioconductor packages require `BiocManager::install()`:

**CRAN:**
* `TwoSampleMR`
* `ieugwasr` (used for OpenGWAS API queries in drug target and downstream analyses)
* `cause`
* `MRPRESSO`
* `tidyverse`
* `data.table`
* `cowplot`
* `ggrepel`
* `knitr`
* `kableExtra`
* `forcats`

**Bioconductor:**
* `GenomicRanges`
* `SNPlocs.Hsapiens.dbSNP144.GRCh37`

## Downloading and Parsing Chromosomes

Most of this is contained in the script `reference-build.sh`.  This was run on GreatLakes.  Key steps include:

* Create Output Directories.  Sets up folders to organize downloaded VCF files and resulting PLINK files.
* Download the 1000 Genomes Sample Panel File.  Retrieves the panel file listing all samples and their population ancestry.
* Extract European (EUR) Sample IDs.  Uses awk to create a list of sample IDs belonging to the EUR super-population.
* Loop Over Chromosomes 1–22. For each autosomal chromosome:  
  * Download Chromosome VCF and Index.  Fetches the compressed VCF file and its index for the chromosome.
  * Subset VCF to EUR Samples. Uses bcftools to keep only EUR individuals in the VCF.
  * Index the Subsetted VCF.  Indexes the new EUR-specific VCF file for efficient access.
  * Subset to only biallelelc SNPs.
  * Convert to PLINK Format.  Uses plink2 to convert the EUR VCF to PLINK binary format (.bed/.bim/.fam).

These temporary files are all stored on GreatLakes or DataDen and not part of this Repository

As a result I get set of PLINK-format genotype files for each chromosome, containing only European ancestry samples from the 1000 Genomes Project—ready for use as an LD reference panel in downstream analyses.

## Summary Statistics Files

### UK Biobank Phenotype Summary Statistics

The manifest for these files including links was found on the [Pan-UK Biobank phenotype manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit?gid=1450719288#gid=1450719288):

* Total Cholesterol (trait 30690; 420607 cases) from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30690-both_sexes-irnt.tsv.bgz
* LDL Cholesterol (trait 30780; 419831 cases) from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30780-both_sexes-irnt.tsv.bgz
* Serum Calcium (trait 30680; 385066 cases) from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30680-both_sexes-irnt.tsv.bgz

These were extracted via zcat for use in the clumping analsyes.

Did LD clumping (using `ld-clumping.slurm` script).  This script has two stages

1. It creates a merged genotype file, after fixing the headers (1000G_EUR_merged) in plink2 format.
2. This reference genome was used to clump the summary statistics for Total Cholesterol, LDL-C, and calcium after some trimming of the file.
3. The ld clumping was done using plink 10000 kb windows with cutoffs at 5E-8 to pick a SNP and 1E-6 as a secondary cutoff.  The clump R2 was set to 0.01.
4. Calculated MAFs for both clumped and all SNPs.

## Filtering and Analysis of Instrument SNPS

SNPS were clumped using plink as above and analyzed using the R script `ld_clump_analysis.qmd`.  We filtered to only include SNPS with a European MAF >= 0.01 and calculated summary statistics for both calcium and cholesterol SNPs.  Wrote out the list of filtered SNPs.

## MR Analyses

### Calcium-Calcium MR Analyses

As a a positive control we analysed UKBB predicted calcium versus MGI-BioVU calcium in `mr-ldlc-calcium.qmd`.  This was a proof of concept that MR analyses could predict calcium levels across these datasets.  This will not be presented in a paper

### Calcium-Cholesterol/LDLC MR Analyses

We next peformed a series of MR and sensitivity analyses testing the directionality of the calcium/cholesterol associations:

- `mr-tc-calcium.qmd` testing the effects of total cholesterol (UKBB) on calcium (MGI-BioVU)
- `mr-ldlc-calcium.qmd` testing the effects of LDL cholesterol (UKBB) on calcium (MGI-BioVU)
- `mr-calcium-tc.qmd` testing the effects of serum calcium (UKBB) on total cholesterol (MGI-BioVU)
- `mr-calcium-ldlc.qmd` testing the effects of serum calcium (UKBB) on LDL-C (MGI-BioVU)

### Drug Target MR Analyses

Our first approach was to use UKBB-based cholesterol SNPs and MGI-BioVU calcium results to evaluate the specific effects of HMGCR, PCKS9, and NPC1L1 on serum calcium but there was very low coverage of the latter two in MGI-BioVU GWAS (see `drug_target_mr_analysis.qmd`), so we repeated this analysis using UKBB-based calcium SNPs (see `drug_target_mr_ukb.qmd`).  The LDL-C exposures were from ebi-a-GCST90025953, a GLGC meta-analysis

### Calcium outcome (completed)
- **Script**: `drug_target_mr_ukb.qmd`
- **Outcome**: UK Biobank serum calcium (Barton 2021, n = 400,792)
- **Finding**: HMGCR cis-instruments significantly affect serum calcium
  (IVW-RE β = 0.18, p < 0.001); PCSK9 informatively null (β = 0.004, p = 0.87,
  F > 2000); NPC1L1 imprecise
- **Interpretation**: Cholesterol → calcium effect is mediated by the
  mevalonate pathway, not LDL-C per se

### BMD outcomes (completed 2026-05-11)
- **Script**: `drug_target_mr_bmd.qmd`
- **Outcomes**: Heel BMD (Morris 2019, n = 426,824), Femoral neck BMD
  (Zheng 2015), Fractures (Dönertaş 2021), Vitamin D (MGI-BioVU; negative
  control)
- **Finding**: HMGCR cis-instruments significantly reduce Heel BMD
  (IVW-RE β = −0.11, p < 0.001 for both LDL-C and total cholesterol);
  PCSK9 null (β ≈ −0.02, p > 0.2); NPC1L1 imprecise
- **Interpretation**: Parallels the calcium finding — the cholesterol → BMD
  pathway is mevalonate-pathway-mediated; LDL-C reduction via the
  PCSK9-LDLR axis does not affect BMD
- **Output**: `MR Results - Drug Target BMD.csv`

### Joint significance
Together these two drug-target MR analyses establish that the
cholesterol → BMD → calcium pathway is HMGCR/mevalonate-mediated at every
step. This is the mechanistic basis for the proposed statin/ezetimibe and
cell-specific Hmgcr knockout experiments.

### Datasets Summary

| Experiment | Exposure/Outcome | Variable | Dataset | ID | n |
|---|---|---|---|---|---|
| mr-tc-calcium | Exposure | Total Cholesterol | UK Biobank | 30690 | 420,607 |
| mr-tc-calcium | Outcome | Serum Calcium | MGI-BioVU LabWAS | | 46,100 |
| mr-ldlc-calcium | Exposure | LDL Cholesterol | UK Biobank | 30780 | 419,831 |
| mr-ldlc-calcium | Outcome | Serum Calcium | MGI-BioVU LabWAS | | 46,100 |
| mr-calcium-tc | Exposure | Serum Calcium | UK Biobank | 30680 | 385,066 |
| mr-calcium-tc | Outcome | Total Cholesterol | MGI-BioVU LabWAS | | 46,100 |
| mr-calcium-ldlc | Exposure | Serum Calcium | UK Biobank | 30680 | 385,066 |
| mr-calcium-ldlc | Outcome | LDL Cholesterol | MGI-BioVU LabWAS | | 46,100 |
| drug_target_mr_analysis | Exposure | Total Cholesterol | GLGC meta-analysis | ebi-a-GCST90025953 | 1,320,016 |
| drug_target_mr_analysis | Outcome | Serum Calcium | MGI-BioVU LabWAS | | 46,100 |
| drug_target_mr_ukb | Exposure | Total Cholesterol | GLGC meta-analysis | ebi-a-GCST90025953 | 1,320,016 |
| drug_target_mr_ukb | Outcome | Serum Calcium | UKB Barton et al. 2021 | ebi-a-GCST90025990 | 400,792 |
| mr-downstream-analyses | Exposure | Total Cholesterol | UK Biobank | 30690 | 420,607 |
| mr-downstream-analyses | Outcome | Vitamin D | MGI-BioVU LabWAS | | 12,250 |
| mr-downstream-analyses | Outcome | Heel BMD | Morris 2019, UK Biobank | ebi-a-GCST006979 | 426,824 |
| mr-downstream-analyses | Outcome | Femoral Neck BMD | Zheng 2015, GEFOS | ieu-a-980 | 32,735 |
| mr-downstream-analyses | Outcome | Fractures | Dönertaş 2021, UK Biobank | ebi-a-GCST90038703 | 484,598 |

### Outcome GWAS

To evaluate the effects of cholesterol on vitamin D and bone mineral density we performed another set of MR analyses using `r mr-downstream-analyses.qmd`.

### Summary Figures and Tables

Each script exported several files including instruments used, instrument summary statistics and a variety of model outputs.  These were combined for making figures and tables in the file `r summary-tables.qmd`.