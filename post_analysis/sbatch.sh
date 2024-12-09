#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=test
#SBATCH --output=out.out
#SBATCH --error=error.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

# Load Conda module (if needed) and activate the environment
source /private/home/rsmerigl/anaconda3/bin/activate r_env


# Esegui lo script Python
python /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_analysis/ancestry_keep_files_creation.py

