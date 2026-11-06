#!/bin/sh
# GEMMA LMM for serum calcium phenotypes — Svenson DO HFD v12
# Run after extract_calcium_phenotype.R has written Calcium2_all.txt etc.
#
# Each phenotype has a different number of non-NA individuals (calcium2=763,
# calcium1=636, non-fasting=409), so GEMMA must build its own kinship and
# eigen decomposition for each — the cholesterol eigen cannot be reused.

#SBATCH --job-name=gemma-calcium
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20g
#SBATCH --time=120:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/logs/%x-%j.log

# ── Primary: 19-week calcium (rank-normalised) ───────────────────────────────
gemma -g Genotypes_all.bimbam -p Calcium2_all.txt -c Covariates.tab \
      -gk 1 -o calcium2_all
gemma -g Genotypes_all.bimbam -p Calcium2_all.txt -c Covariates.tab \
      -k output/calcium2_all.cXX.txt -eigen -o calcium2_all
gemma -g Genotypes_all.bimbam -p Calcium2_all.txt -c Covariates.tab \
      -d output/calcium2_all.eigenD.txt -u output/calcium2_all.eigenU.txt \
      -lmm 1 -o calcium2_all

# ── Secondary: 10-week calcium ───────────────────────────────────────────────
gemma -g Genotypes_all.bimbam -p Calcium1_all.txt -c Covariates.tab \
      -gk 1 -o calcium1_all
gemma -g Genotypes_all.bimbam -p Calcium1_all.txt -c Covariates.tab \
      -k output/calcium1_all.cXX.txt -eigen -o calcium1_all
gemma -g Genotypes_all.bimbam -p Calcium1_all.txt -c Covariates.tab \
      -d output/calcium1_all.eigenD.txt -u output/calcium1_all.eigenU.txt \
      -lmm 1 -o calcium1_all

# ── Sensitivity: non-fasting calcium ─────────────────────────────────────────
gemma -g Genotypes_all.bimbam -p Calcium_nonfasting_all.txt -c Covariates.tab \
      -gk 1 -o calcium_nonfasting_all
gemma -g Genotypes_all.bimbam -p Calcium_nonfasting_all.txt -c Covariates.tab \
      -k output/calcium_nonfasting_all.cXX.txt -eigen -o calcium_nonfasting_all
gemma -g Genotypes_all.bimbam -p Calcium_nonfasting_all.txt -c Covariates.tab \
      -d output/calcium_nonfasting_all.eigenD.txt \
      -u output/calcium_nonfasting_all.eigenU.txt \
      -lmm 1 -o calcium_nonfasting_all
