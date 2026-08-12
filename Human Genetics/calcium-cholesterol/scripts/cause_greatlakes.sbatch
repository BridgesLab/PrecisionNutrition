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
#   ~/logs/                              MUST EXIST BEFORE SUBMITTING — slurm opens
#                                        --output before this script runs
#
# Submit both arms (they are independent, so run them in parallel):
#   mkdir -p ~/logs
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
# %u = username, %x = job name, %j = job id. These paths are resolved by slurm
# BEFORE the script runs, so the directory must already exist:  mkdir -p ~/logs
#SBATCH --output=/home/%u/logs/%x_%j.out
#SBATCH --error=/home/%u/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davebrid@umich.edu

set -euo pipefail

ARM="${ARM:-overlapping}"
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR}}"

# Module names drift between cluster upgrades, and `set -e` would turn a failed
# `module load` into a silent job death. Try, and carry on if unavailable.
module purge 2>/dev/null || true
for m in R gcc; do
  module load "$m" 2>/dev/null && echo "loaded module: $m" || echo "no module: $m"
done

# Great Lakes uses a hierarchical module tree: the Bioinformatics collection has
# to be loaded before plink is visible.
module load Bioinformatics 2>/dev/null && echo "loaded module: Bioinformatics" \
  || echo "no module: Bioinformatics"

# plink/1.9 SPECIFICALLY. ieugwasr::ld_clump() uses plink 1.9's --clump syntax
# and output format; plink 2.x changed both, and would fail at the pruning step
# after the expensive parameter estimation has already run.
module load plink/1.9 2>/dev/null && echo "loaded module: plink/1.9" \
  || echo "WARNING: could not load plink/1.9"

command -v Rscript >/dev/null || { echo "FATAL: no Rscript after module load. Try 'module spider R'."; exit 1; }

PLINK_BIN="$(command -v plink || true)"
export PLINK_BIN
if [[ -z "${PLINK_BIN}" ]]; then
  echo "no plink on PATH — resolve_plink() will try genetics.binaRies"
else
  PLINK_VER="$(plink --version 2>/dev/null | head -1)"
  echo "plink version : ${PLINK_VER}"
  case "${PLINK_VER}" in
    *v1.9*|*v1.90*) : ;;
    *) echo "FATAL: ld_clump() requires plink 1.9, got: ${PLINK_VER}"
       echo "       module load Bioinformatics && module load plink/1.9"
       exit 1 ;;
  esac
fi

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
