#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=BMI_SAS
#SBATCH --output=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/BMI/output_BMI_SAS.txt
#SBATCH --error=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/BMI/error_BMI_SAS.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/admixture
python /private/home/rsmerigl/codes/Admixture_mapping/association_execution/BMI/association_execution_BMI_SAS.py