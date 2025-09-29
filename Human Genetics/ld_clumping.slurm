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

### GET AND PREPARE SUMMARY STATISTICS

# LDL-C
CAL_FILE="biomarkers-30680-both_sexes-irnt.tsv.bgz"

# Download if not exists
if [ ! -f "$CAL_FILE" ]; then
    wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/$CAL_FILE
fi

# Extract variant ID, beta, SE, and P-value if uncompressed file missing
if [ ! -f biomarkers-30680-both_sexes-irnt.tsv ]; then
    zcat "$CAL_FILE" | awk '
    BEGIN { OFS="\t" }
    NR==1 { print "VARIANT", "BETA", "SE", "P" }
    NR>1 {
        p = 10^(-$8);
        variant = $1 ":" $2 ":" $3 ":" $4;
        print variant, $6, $7, p
    }' > biomarkers-30680-both_sexes-irnt.tsv
fi

CAL_FILE="biomarkers-30780-both_sexes-irnt.tsv.bgz"

# Download if not exists
if [ ! -f "$CAL_FILE" ]; then
    wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/$CAL_FILE
fi

# Extract variant ID, beta, SE, and P-value if uncompressed file missing
if [ ! -f biomarkers-30780-both_sexes-irnt.tsv ]; then
    zcat "$CAL_FILE" | awk '
    BEGIN { OFS="\t" }
    NR==1 { print "VARIANT", "BETA", "SE", "P" }
    NR>1 {
        p = 10^(-$8);
        variant = $1 ":" $2 ":" $3 ":" $4;
        print variant, $6, $7, p
    }' > biomarkers-30780-both_sexes-irnt.tsv
fi

CAL_FILE="biomarkers-30690-both_sexes-irnt.tsv.bgz"

# Download if not exists
if [ ! -f "$CAL_FILE" ]; then
    wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/$CAL_FILE
fi

# Extract variant ID, beta, SE, and P-value if uncompressed file missing
if [ ! -f biomarkers-30690-both_sexes-irnt.tsv ]; then
    zcat "$CAL_FILE" | awk '
    BEGIN { OFS="\t" }
    NR==1 { print "VARIANT", "BETA", "SE", "P" }
    NR>1 {
        p = 10^(-$8);
        variant = $1 ":" $2 ":" $3 ":" $4;
        print variant, $6, $7, p
    }' > biomarkers-30690-both_sexes-irnt.tsv
fi

# LD CLUMPING

# Total cholesterol
plink2 \
  --pfile 1000G_EUR_merged \
  --clump biomarkers-30690-both_sexes-irnt.tsv \
  --clump-snp-field VARIANT \
  --clump-field P \
  --clump-p1 5e-8 \
  --clump-p2 1e-6 \
  --clump-r2 0.01 \
  --clump-kb 10000 \
  --threads 8 \
  --out biomarkers-30690-both_sexes-irnt-results

# LDL cholesterol
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
awk 'NR > 1 {print $3}' biomarkers-30690-both_sexes-irnt-results.clumps > biomarkers-30690-both_sexes-irnt-results-clumped-snps.txt


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

   plink2 \
  --pfile 1000G_EUR_merged\
  --extract biomarkers-30690-both_sexes-irnt-results-clumped-snps.txt \
  --freq \
  --threads 8 \
  --out biomarkers-30690-both_sexes-irnt-results-clumped-snps

# For all-snps for 30680
awk 'NR>1 {print $3; gsub(/,/, "\n", $14); print $14}' biomarkers-30680-both_sexes-irnt-results.clumps \
  | sed '/^\s*$/d' \
  | sort -u \
  > biomarkers-30680-all-snps.txt

# For all snps for 30780
awk 'NR>1 {print $3; gsub(/,/, "\n", $14); print $14}' biomarkers-30780-both_sexes-irnt-results.clumps \
  | sed '/^\s*$/d' \
  | sort -u \
  > biomarkers-30780-all-snps.txt

  # For all snps for 30690
awk 'NR>1 {print $3; gsub(/,/, "\n", $14); print $14}' biomarkers-30690-both_sexes-irnt-results.clumps \
  | sed '/^\s*$/d' \
  | sort -u \
  > biomarkers-30690-all-snps.txt

#calculated MAF for all SNPs
plink2 \
  --pfile 1000G_EUR_merged \
  --extract biomarkers-30680-all-snps.txt \
  --freq \
  --threads 8 \
  --out biomarkers-30680-all-snps

plink2 \
  --pfile 1000G_EUR_merged \
  --extract biomarkers-30780-all-snps.txt \
  --freq \
  --threads 8 \
  --out biomarkers-30780-all-snps

  plink2 \
  --pfile 1000G_EUR_merged \
  --extract biomarkers-30690-all-snps.txt \
  --freq \
  --threads 8 \
  --out biomarkers-30690-all-snps