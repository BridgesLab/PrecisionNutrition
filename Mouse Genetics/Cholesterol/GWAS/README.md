# GWAS

Created bimbam genotype files and annotation files (Annotation.csv, Covariates.[diet/sex/gxe].tab, and Cholesterol_[all/ncd/hfd]) using `bimbam-generation.sh` from the version 12 dataset downloaded `core.Svenson_DO_HFD.v12.Rdata ` on November 3, 2023.
Ran GEMMA for several analyses or cholesterol analyses using `gemma-all.sh`

In order to generate PRS and clumping had to re-generate genotype files for Plink format (bed/bim/fam files).  This was done using **Plink-generation.qmd** within the prior script which had to ignore the dosage data and just round each snp up to 0/1/2 (Plink 1.9 does not use dosage data, Plink 2+ didnt have ld calculations set up yet.

The GWAS analyses were run for these groups:

* all - adjusting for diet and sex
* ncd - only NCD, adjusting for sex
* hfd - only HFD, adjusting for sex
* males - only males, adjusting for diet
* females - only females, adjusting for diet
* gxe - adjusting for sex looking for gxe effects (gene x diet), this isnt working.

To run bslmm used `gemma-bslmm.sh` using sex and diet adjusted covariates
* bslmm - adjusting cholesterol for diet and sex first, then doing bslmm analyses

From the gamma parameter (posterior inclusion probability) from the cholesterol_bslmm.param.txt it appears that the highest probability is only 0.0098 (chromosome 4), so this suggests a highly polygenic trait.

These QTLs were analysed using the `qtl-analysis.Rmd` script.

# TWAS

Used the downloaded RDS files from the Svensson DO HFD Viewer to run TWAS.  Used DESeq2 to analyze the raw counts table using the script `twas.sh` which ran the R script **liver-mRNA-DESeq.qmd**

# eQTL Identification

Used gemma to analyse transcript-level eQTLs for selected genes (see eQTL folder).

This used the script `eQTL_additive.sh` and was an additive model (sex + diet).  These results were analysed using the R script `eqtl-analysis.Rmd`

# Linkage disequillibrium estimates

First converted the GEMMA assoc.txt files into the proper format for plink.  This used the `assoc_conversion.py` script to generate .assoc files in the proper format:

python assoc_conversion.py ../output/cholesterol_hfd.assoc.txt 
python assoc_conversion.py ../output/cholesterol_ncd.assoc.txt 
python assoc_conversion.py ../output/cholesterol_females.assoc.txt 
python assoc_conversion.py ../output/cholesterol_males.assoc.txt 
python assoc_conversion.py ../output/cholesterol_all.assoc.txt 

## Calculating clumps and LD tables

Clumping and ld calculations were done with the `plink-postprocessing.sh` script

Clumps were annotated using `clumping-summary.qmd` and analysed using `clumps-twas-integration-diet.qmd` and `clumps-twas-integration-sex.qmd`.

## To Visualize the Nominated QTLs as Barplots

The `snps-of-interest.Rmd` script generates barplots of nominated SNPs.  These SNPs are pulled from the main genotype file from the main genotype file by running the command, which uses a list of SNPs (in `SNPs_of_interest.txt `) to make a shorter genotypes file.  These are pulled from lead snps from the clumps.

