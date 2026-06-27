#!/bin/sh
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=bimbam-generation
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=140g
#SBATCH --time=40:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/logs/%x-%j.log


module load R

quarto render Bimbam-generation.qmd
quarto render SNP-Annotation.qmd
quarto render Phenotype-setup.qmd
