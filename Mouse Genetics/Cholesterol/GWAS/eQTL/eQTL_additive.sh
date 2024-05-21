#!/bin/sh
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=eQTL_additive
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=15g
#SBATCH --time=40:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/logs/%x-%j.log
#SBATCH --error=/home/%u/error_logs/%x-%j.err

GENES=('ENSMUSG00000005681' 'ENSMUSG00000037936')
GENOTYPES='../Genotypes_all.bimbam'
COVARIATES='../Covariates.tab'
MRNA='../mRNA_Expression_all.txt'

module load python
pip install csvkit

for gene in "${GENES[@]}";
do
./select_gene.sh '../mRNA_Expression_all.txt' ${gene} output/${gene}.txt
gemma -g ../Genotypes_all.bimbam -p output/${gene}.txt -c ../Covariates.tab -gk 1 -o ${gene}
gemma -g ../Genotypes_all.bimbam -p output/${gene}.txt -c ../Covariates.tab -k output/${gene}.cXX.txt -eigen -o ${gene}
gemma -g ../Genotypes_all.bimbam -p output/${gene}.txt -c ../Covariates.tab -d output/${gene}.eigenD.txt -u output/${gene}.eigenU.txt -lmm 1 -o ${gene}
rm output output/${gene}.cXX.txt
rm output/${gene}.eigenD.txt
rm output/${gene}.eigenU.txt

done
