#!/bin/bash
#SBATCH --job-name calcium_ld_clumping
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G 
#SBATCH --time=1:00:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/%u/logs/%x-%j.log
#SBATCH --error=/home/%u/error_logs/error-%x-%j.log

module load Bioinformatics
module load samtools
module load bcftools
module load plink2

echo "Running from $(pwd)"

# LD CLUMPING

# LDL
plink2 \
  --pfile 1000G_EUR_merged \
  --clump biomarkers-30780-both_sexes-irnt.tsv \
  --clump-snp-field VARIANT \
  --clump-field P \
  --clump-p1 5e-8 \
  --clump-p2 1e-6 \
  --clump-r2 0.01 \
  --clump-kb 10000 \
  --threads 8 \
  --out biomarkers-30780-both_sexes-irnt-results

# Calcium
plink2 \
  --pfile 1000G_EUR_merged \
  --clump biomarkers-30680-both_sexes-irnt.tsv \
  --clump-snp-field VARIANT \
  --clump-field P \
  --clump-p1 5e-8 \
  --clump-p2 1e-6 \
  --clump-r2 0.01 \
  --clump-kb 10000 \
  --threads 8 \
  --out biomarkers-30680-both_sexes-irnt-results

# created list of clumped SNPs
awk 'NR > 1 {print $3}' biomarkers-30680-both_sexes-irnt-results.clumps > biomarkers-30680-both_sexes-irnt-results-clumped-snps.txt
awk 'NR > 1 {print $3}' biomarkers-30780-both_sexes-irnt-results.clumps > biomarkers-30780-both_sexes-irnt-results-clumped-snps.txt

# calculated MAF for clumped snps
plink2 \
  --pfile 1000G_EUR_merged \
  --extract biomarkers-30680-both_sexes-irnt-results-clumped-snps.txt \
  --freq \
  --threads 8 \
  --out biomarkers-30680-both_sexes-irnt-results-clumped-snps

 plink2 \
  --pfile 1000G_EUR_merged\
  --extract biomarkers-30780-both_sexes-irnt-results-clumped-snps.txt \
  --freq \
  --threads 8 \
  --out biomarkers-30780-both_sexes-irnt-results-clumped-snps
