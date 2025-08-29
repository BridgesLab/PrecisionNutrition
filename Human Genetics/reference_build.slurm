#!  /bin/bash
#SBATCH --job-name calcium_reference_build
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G 
#SBATCH --time=4:00:00
#SBATCH --account=davebrid0
#SBATCH --partition=standard
#SBATCH --mail-user=davebrid@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/%u/logs/%x-%j.log
#SBATCH --error=/home/%u/error_logs/error-%x-%j.log

module load Bioinformatics
module load samtools

echo "Running from $(pwd)"

./reference_build.sh
