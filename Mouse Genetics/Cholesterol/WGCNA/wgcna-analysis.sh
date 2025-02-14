#!/bin/sh
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=wgcna
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=7g
#SBATCH --time=90:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/logs/%x-%j.log
#SBATCH --error=/home/%u/error_logs/%x-%j.log

module load R

quarto render Liver-WGCNA-Construction.Qmd 
quarto render Liver-WGCNA-Cholesterol.Qmd