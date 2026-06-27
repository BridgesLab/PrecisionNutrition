#!/bin/sh
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=gemma-all
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
 
gemma -g Genotypes_all.bimbam -p Cholesterol_all.txt -c Covariates.tab -gk 1 -o cholesterol_all
gemma -g Genotypes_all.bimbam -p Cholesterol_all.txt -c Covariates.tab -k output/cholesterol_all.cXX.txt -eigen -o cholesterol_all
gemma -g Genotypes_all.bimbam -p Cholesterol_all.txt -c Covariates.tab -d output/cholesterol_all.eigenD.txt -u output/cholesterol_all.eigenU.txt -lmm 1 -o cholesterol_all


gemma -g Genotypes_all.bimbam -p Cholesterol_HFD.txt -c Covariates.sex.tab -gk 1 -o cholesterol_hfd
gemma -g Genotypes_all.bimbam -p Cholesterol_HFD.txt -c Covariates.sex.tab -k output/cholesterol_hfd.cXX.txt -eigen -o cholesterol_hfd
gemma -g Genotypes_all.bimbam -p Cholesterol_HFD.txt -c Covariates.sex.tab -d output/cholesterol_hfd.eigenD.txt -u output/cholesterol_hfd.eigenU.txt -lmm 1 -o cholesterol_hfd


gemma -g Genotypes_all.bimbam -p Cholesterol_NCD.txt -c Covariates.sex.tab -gk 1 -o cholesterol_ncd
gemma -g Genotypes_all.bimbam -p Cholesterol_NCD.txt -c Covariates.sex.tab -k output/cholesterol_ncd.cXX.txt -eigen -o cholesterol_ncd
gemma -g Genotypes_all.bimbam -p Cholesterol_NCD.txt -c Covariates.sex.tab -d output/cholesterol_ncd.eigenD.txt -u output/cholesterol_ncd.eigenU.txt -lmm 1 -o cholesterol_ncd
