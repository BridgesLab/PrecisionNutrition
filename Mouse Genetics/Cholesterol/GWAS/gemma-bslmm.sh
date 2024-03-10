#!/bin/sh
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=gemma-bslmm
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20g
#SBATCH --time=100:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/logs/%x-%j.log
 
gemma -g Genotypes_all.bimbam -p Cholesterol_adj.txt -bslmm 1 -o cholesterol_bslmm
