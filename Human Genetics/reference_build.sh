#!/bin/bash

set -e

# Create directories for organization
mkdir -p 1000G_VCF 1000G_EUR_PLINK

# Download the panel file from EBI, fallback to UW Beagle mirror if needed
PANEL_EBI="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
PANEL_UW="http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel"

echo "Downloading 1000 Genomes sample panel file..."
if [ ! -s 1000G_VCF/integrated_call_samples_v3.20130502.ALL.panel ]; then
    curl --max-time 300 --retry 3 --retry-delay 5 -f -o 1000G_VCF/integrated_call_samples_v3.20130502.ALL.panel "$PANEL_EBI" || {
        echo "EBI download failed, trying UW Beagle mirror..."
        curl --max-time 300 --retry 3 --retry-delay 5 -f -o 1000G_VCF/integrated_call_samples_v3.20130502.ALL.panel "$PANEL_UW" || {
            echo "Failed to download panel file from both sources. Exiting."
            exit 1
        }
    }
fi

if [ ! -s 1000G_VCF/integrated_call_samples_v3.20130502.ALL.panel ]; then
    echo "Panel file is missing or empty. Exiting."
    exit 1
fi

# Extract EUR sample IDs (POSIX awk)
awk '$3 == "EUR" {print $1}' 1000G_VCF/integrated_call_samples_v3.20130502.ALL.panel > 1000G_VCF/eur_samples.txt

# Function to check disk space (at least 10G free)
check_disk_space() {
    FREE=$(df -Pk . | awk 'NR==2 {print $4}')
    # 10G = 10*1024*1024 KB = 10485760 KB
    if [ "$FREE" -lt 10485760 ]; then
        echo "Warning: Less than 10GB free disk space. Aborting."
        exit 1
    fi
}

# Loop over chromosomes 1-22
for chr in $(seq 1 22); do
    echo "Processing chromosome $chr..."

    VCF_FN="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    TBI_FN="${VCF_FN}.tbi"
    VCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${VCF_FN}"
    TBI_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${TBI_FN}"

    # Skip if already processed
    if [ -s 1000G_EUR_PLINK/chr${chr}_EUR.bed ]; then
        echo "PLINK files for chr${chr} already exist, skipping."
        continue
    fi

    # Check disk space before download
    check_disk_space

    # Download VCF if not already present
    if [ ! -s 1000G_VCF/${VCF_FN} ]; then
        echo "Downloading ${VCF_FN}..."
        if ! curl --max-time 1800 --retry 5 --retry-delay 10 -f -o 1000G_VCF/${VCF_FN} "$VCF_URL"; then
            echo "VCF for chr${chr} not found or download failed, skipping..."
            continue
        fi
    else
        echo "${VCF_FN} already downloaded."
    fi

    # Download TBI if not already present
    if [ ! -s 1000G_VCF/${TBI_FN} ]; then
        echo "Downloading ${TBI_FN}..."
        if ! curl --max-time 600 --retry 5 --retry-delay 10 -f -o 1000G_VCF/${TBI_FN} "$TBI_URL"; then
            echo "TBI for chr${chr} not found or download failed, skipping..."
            continue
        fi
    else
        echo "${TBI_FN} already downloaded."
    fi

    # Subset to EUR samples if not already done
    if [ ! -s 1000G_VCF/chr${chr}.EUR.vcf.gz ]; then
        echo "Subsetting chr${chr} to EUR samples..."
        bcftools view -S 1000G_VCF/eur_samples.txt -Oz \
            -o 1000G_VCF/chr${chr}.EUR.vcf.gz 1000G_VCF/${VCF_FN}
        tabix -p vcf 1000G_VCF/chr${chr}.EUR.vcf.gz
    else
        echo "chr${chr}.EUR.vcf.gz already exists."
    fi

    # Filter to biallelic SNPs only before PLINK conversion
    if [ ! -s 1000G_VCF/chr${chr}.EUR.biallelic.vcf.gz ]; then
        echo "Filtering chr${chr}.EUR.vcf.gz to biallelic SNPs only..."
        bcftools view -m2 -M2 -v snps -Oz \
            -o 1000G_VCF/chr${chr}.EUR.biallelic.vcf.gz 1000G_VCF/chr${chr}.EUR.vcf.gz
        tabix -p vcf 1000G_VCF/chr${chr}.EUR.biallelic.vcf.gz
    else
        echo "chr${chr}.EUR.biallelic.vcf.gz already exists."
    fi

done

echo "All chromosomes processed! Your European LD reference panel is ready in 1000G_EUR_PLINK/"
