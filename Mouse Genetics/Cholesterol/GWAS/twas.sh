#!/bin/sh
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=do-hfd-twas
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=1g
#SBATCH --time=65:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/logs/%x-%j.log
#SBATCH --error=/home/%u/error_logs/%x-%j.log

module load R

#quarto render liver-mRNA-DESeq.qmd
quarto render liver-mRNA-lm-analyses.qmd

