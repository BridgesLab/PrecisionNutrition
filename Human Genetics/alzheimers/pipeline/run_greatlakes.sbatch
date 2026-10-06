#!/bin/bash
# run_greatlakes.sbatch — run the targets pipeline on Great Lakes.
#
# This job is the CONTROLLER: it runs targets::tar_make(), which does light work itself and
# submits heavy targets (MRAID, CAUSE, the SlopeHunter checks) as separate Slurm jobs through
# crew.cluster, using the resources in pipeline/config.yml (greatlakes:). Submit from alzheimers/:
#   mkdir -p ~/logs logs && sbatch pipeline/run_greatlakes.sbatch
# Re-submitting resumes: targets that are already up to date are skipped.

#SBATCH --job-name=ad_robust_mr
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/home/%u/logs/%x_%j.out
#SBATCH --error=/home/%u/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davebrid@umich.edu

set -uo pipefail

# Module names drift between cluster upgrades; try each and report rather than die silently.
module purge 2>/dev/null || true
for m in R gcc; do
  module load "$m" 2>/dev/null && echo "loaded module: $m" || echo "no module: $m (try: module spider $m)"
done

# plink2 from ~/bin, NOT a module: Great Lakes' newest plink2 module (2.00a4, Jan 2023) has no
# --clump, which CAUSE's LD pruning needs. Install: see pipeline/README.md. Fail now, not after
# an hour of CAUSE parameter estimation, if the binary is missing or cannot clump.
PLINK2_BIN="${PLINK2_BIN:-$HOME/bin/plink2}"
export PLINK2_BIN
[[ -x "${PLINK2_BIN}" ]] || { echo "FATAL: ${PLINK2_BIN} not found or not executable"; exit 1; }
"${PLINK2_BIN}" --help clump 2>&1 | grep -q -- "--clump-p1" \
  || { echo "FATAL: ${PLINK2_BIN} has no --clump ($("${PLINK2_BIN}" --version | head -1))"; exit 1; }
export PIPELINE_PROFILE=greatlakes

cd "${SLURM_SUBMIT_DIR}"
echo "Started $(date) on $(hostname); plink2 = ${PLINK2_BIN}"
# Notebook targets (nb_*) are excluded: they pull from the OpenGWAS API, which needs a token on the
# cluster and was down when this was written. Render them locally instead.
Rscript -e 'targets::tar_make(names = !tidyselect::starts_with("nb_"))'
status=$?
Rscript -e 'print(targets::tar_progress_summary()); m <- targets::tar_meta(fields = error, complete_only = TRUE); print(m[!is.na(m$error), ])'
echo "Finished $(date), tar_make exit ${status}"
exit ${status}
