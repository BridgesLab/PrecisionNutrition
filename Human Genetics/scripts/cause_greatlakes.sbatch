#!/bin/bash
# cause_greatlakes.sbatch -------------------------------------------------------
# CAUSE is the slowest piece of this arm: est_cause_params() on ~1M variants plus
# the three-model fit runs ~1-2 h at 16-32 GB on the HapMap3-restricted input.
# That will also run on a laptop overnight; the cluster just makes it painless and
# lets both arms go in parallel. MR-APSS takes minutes and is NOT sent here.
#
# WHAT MUST BE ON THE CLUSTER (this file alone is not enough):
#   scripts/cause_greatlakes.sbatch      this wrapper
#   scripts/run_cause_greatlakes.R       the actual job
#   R/robust_mr_helpers.R                sourced by the above
#   config_robust_mr.yml                 paths and thresholds
#   raw_data/robust_mr/cache/*_cause.rds harmonised inputs from robust_mr_prep.qmd
#   raw_data/reference/1kg_eur/EUR.{bed,bim,fam}   LD panel for pruning
#   logs/                                MUST EXIST BEFORE SUBMITTING — slurm opens
#                                        --output before this script runs
#
# Submit both arms (they are independent, so run them in parallel):
#   mkdir -p logs
#   sbatch --export=ALL,ARM=overlapping  scripts/cause_greatlakes.sbatch
#   sbatch --export=ALL,ARM=overlap-free scripts/cause_greatlakes.sbatch
#
# Then rsync the cache back and knit robust_mr_cause.qmd locally — it picks up
# cause_fit_<arm>.rds automatically and skips the expensive chunks.

#SBATCH --job-name=cause_ldl_bmd
#SBATCH --account=davebrid0              # your Great Lakes slurm account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
# 48G/6h is sized for the HapMap3-restricted input produced by robust_mr_prep.qmd
# (~1M variants). If you drop the HM3 restriction to run CAUSE's native
# full-genome-wide workflow (~8M variants), raise this to --mem=96G --time=12:00:00.
#SBATCH --output=logs/cause_%x_%j.out
#SBATCH --error=logs/cause_%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davebrid@umich.edu

set -euo pipefail

ARM="${ARM:-overlapping}"
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR}}"

module purge
module load R/4.4.0            # `module spider R` to see what Great Lakes currently offers
module load gcc
module load plink/1.9

# clump_local() picks this up via resolve_plink(); genetics.binaRies is often not
# installed on the cluster, so hand it the module-provided binary explicitly.
export PLINK_BIN="$(command -v plink)"

echo "Host          : $(hostname)"
echo "Arm           : ${ARM}"
echo "Project root  : ${PROJECT_ROOT}"
echo "R             : $(command -v Rscript)"
echo "plink         : ${PLINK_BIN}"
echo "Started       : $(date)"

# Fail fast and legibly if the transfer was incomplete, rather than 40 minutes in.
for f in scripts/run_cause_greatlakes.R R/robust_mr_helpers.R config_robust_mr.yml; do
  [[ -f "${PROJECT_ROOT}/${f}" ]] || { echo "MISSING: ${f}"; exit 1; }
done
for f in glgc2021_ldl willer2013_ldl ebmd; do
  p="${PROJECT_ROOT}/raw_data/robust_mr/cache/${f}_cause.rds"
  [[ -f "${p}" ]] || echo "WARNING: missing ${p} (needed only if ARM uses it)"
done

# One-time package setup (safe to re-run; installs into your personal library):
#   Rscript -e 'install.packages(c("remotes","tidyverse","data.table","yaml","here","ieugwasr"), repos="https://cloud.r-project.org")'
#   Rscript -e 'remotes::install_github(c("stephenslab/mixsqp","stephens999/ashr","jean997/cause"))'

cd "${PROJECT_ROOT}"
Rscript scripts/run_cause_greatlakes.R "${ARM}" "${PROJECT_ROOT}"

echo "Finished      : $(date)"
