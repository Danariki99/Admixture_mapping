#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=test
#SBATCH --output=out_chr%a_%A.out
#SBATCH --error=error_chr%a_%A.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

# Controlla se è stato fornito un argomento
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch submit_job.sh <chromosome_number>"
    exit 1
fi

CHROMOSOME=$1

# Load Conda module (if needed) and activate the environment
source /private/home/rsmerigl/anaconda3/bin/activate r_env

# Esegui lo script Python, passando il cromosoma come argomento
python /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_analysis/ancestry_counts_files.py "$CHROMOSOME" > "output_chr${CHROMOSOME}.log" 2> "error_chr${CHROMOSOME}.log"
