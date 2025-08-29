# Analysis of Summary Statistics for Calcium and LDL-C

We wanted to be able to download summary statistics for some traits, and do clumping analysis using plink2

## Software Requirements

* Downloaded and installed plink2 (v2.00a5.12LM AVX2 AMD (25 Jun 2024))

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

Downloaded from UK Biobank and did LD clumping (using `ld-clumping.slurm` script).  This script has two stages

1. It creates a merged genotype file, after fixing the headers (1000G_EUR)merged) in plnk2 format.
2. This reference genome was used to clump the summary statistics for the LDL-C (biomarkers-30780-both_sexes-irnt.tsv) and calcium (biomarkers-30680-both_sexes-irnt.tsv) GWAS downloaded from UK Biobank.  Again these large tsv files are stored on GreatLakes or DataDen and not in this repository.
2. The ld clumping was done using plink 10000 kb windows with cutoffs at 5E-8 to pick a SNP and 1E-6 as a secondary cutoff.  The clump R2 is set to 0.01.

## Analysis of Instrument SNPS

SNPS were clumped using plink as above and analyzed using the R script `ld_clump_analysis.qmd`