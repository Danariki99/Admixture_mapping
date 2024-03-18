#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=no_BMI
#SBATCH --output=/private/home/rsmerigl/codes/association_execution/output_no_BMI.txt
#SBATCH --error=/private/home/rsmerigl/codes/association_execution/error_no_BMI.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/admixture
python /private/home/rsmerigl/codes/association_execution/association_execution_no_BMI.py