#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=pipeline_execution
#SBATCH --output=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/execution_output.txt
#SBATCH --error=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/execution_error.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/e_env
python /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/pre-processing/pre_processing.py
./ /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/association_execution/job_submission_BMI.sh
python /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post-processing/post_processing.py