#!/bin/bash
#SBATCH --job-name calcium_ld_clumping_setup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G 
#SBATCH --time=3:00:00
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

REF_FASTA="reference/human_g1k_v37.fasta"  # Replace with actual path

# Process each chromosome
for chr in {1..22}; do
  echo "Processing chromosome ${chr}..."

  RAW_VCF="1000G_VCF/chr${chr}.EUR.biallelic.vcf.gz"
  FIXED_VCF="1000G_VCF/chr${chr}.fixed.vcf.gz"
  NORMALIZED_VCF="1000G_VCF/chr${chr}.normalized.vcf.gz"
  FILTERED_VCF="1000G_VCF/chr${chr}.filtered.vcf.gz"
  RENAMED_VCF="1000G_VCF/chr${chr}.renamed.vcf.gz"
  PLINK_OUT="1000G_EUR_PLINK/chr${chr}_EUR"

  # Reheader to fix malformed INFO tags if needed
  bcftools view -h "$RAW_VCF" > header_chr${chr}.txt
  sed -i 's/ID=SVLEN,Number=\.,/ID=SVLEN,Number=A,/' header_chr${chr}.txt
  bcftools reheader -h header_chr${chr}.txt -o "$FIXED_VCF" "$RAW_VCF"
  bcftools index "$FIXED_VCF"
  rm header_chr${chr}.txt

# Normalize with reference genome
bcftools norm \
  -f $REF_FASTA \
  -m-any \
  -Oz -o 1000G_VCF/chr${chr}.normalized.vcf.gz \
  1000G_VCF/chr${chr}.fixed.vcf.gz
bcftools index 1000G_VCF/chr${chr}.normalized.vcf.gz

# Fix REF/ALT inconsistencies using the reference
bcftools +fixref \
  1000G_VCF/chr${chr}.normalized.vcf.gz \
  -- \
  -f $REF_FASTA \
  -Oz -o 1000G_VCF/chr${chr}.fixed_ref.vcf.gz
bcftools index 1000G_VCF/chr${chr}.fixed_ref.vcf.gz

# Add proper ID in format CHR:POS:REF:ALT
bcftools annotate \
  -x ID \
  -I +'%CHROM:%POS:%REF:%ALT' \
  -Oz -o 1000G_VCF/chr${chr}.renamed.vcf.gz \
  1000G_VCF/chr${chr}.fixed_ref.vcf.gz
bcftools index 1000G_VCF/chr${chr}.renamed.vcf.gz

# Convert to PLINK2 format with correct IDs
plink2 \
  --vcf 1000G_VCF/chr${chr}.renamed.vcf.gz \
  --make-pfile \
  --out 1000G_EUR_PLINK/chr${chr}_EUR
done

# Create a list to merge the chromosomes
for chr in {1..22}; do
  echo "1000G_EUR_PLINK/chr${chr}_EUR"
done > pmerge_list.txt

# Merge all chromosomes
plink2 \
  --pmerge-list pmerge_list.txt \
  --make-pfile \
  --threads 4 \
  --out 1000G_EUR_merged