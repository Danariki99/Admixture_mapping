#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=no_BMI_AFR
#SBATCH --output=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/no_BMI/output_no_BMI_AFR.txt
#SBATCH --error=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/no_BMI/error_no_BMI_AFR.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/admixture
python /private/home/rsmerigl/codes/Admixture_mapping/association_execution/no_BMI/association_execution_no_BMI_AFR.py