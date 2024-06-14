#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=pipeline_execution
#SBATCH --output=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/execution_output.txt
#SBATCH --error=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/execution_error.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

# Check if a dataset argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset>"
    exit 1
fi

# Assign the first argument to the variable 'dataset'
dataset=$1

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/r_env

# Pass the dataset variable to the Python scripts
python pre_processing/pre_processing.py $dataset
./association_execution/job_submission_BMI.sh $dataset
python post_processing/post_processing.py $dataset