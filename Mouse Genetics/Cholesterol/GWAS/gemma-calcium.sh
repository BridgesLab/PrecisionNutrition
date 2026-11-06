#!/bin/sh
# GEMMA LMM for serum calcium phenotypes — Svenson DO HFD v12
# Run after extract_calcium_phenotype.R has written Calcium2_all.txt etc.
#
# Reuses the cholesterol_all kinship/eigen decomposition (same genotypes,
# same mice) to avoid recomputing the kinship matrix.

#SBATCH --job-name=gemma-calcium
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20g
#SBATCH --time=40:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/logs/%x-%j.log

# Primary phenotype: 19-week calcium (rank-normalised), reusing chol eigen
gemma -g Genotypes_all.bimbam \
      -p Calcium2_all.txt \
      -c Covariates.tab \
      -d output/cholesterol_all.eigenD.txt \
      -u output/cholesterol_all.eigenU.txt \
      -lmm 1 \
      -o calcium2_all

# Secondary: 10-week calcium
gemma -g Genotypes_all.bimbam \
      -p Calcium1_all.txt \
      -c Covariates.tab \
      -d output/cholesterol_all.eigenD.txt \
      -u output/cholesterol_all.eigenU.txt \
      -lmm 1 \
      -o calcium1_all

# Non-fasting calcium (sensitivity)
gemma -g Genotypes_all.bimbam \
      -p Calcium_nonfasting_all.txt \
      -c Covariates.tab \
      -d output/cholesterol_all.eigenD.txt \
      -u output/cholesterol_all.eigenU.txt \
      -lmm 1 \
      -o calcium_nonfasting_all
