#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name=test
#SBATCH --output=/private/home/rsmerigl/codes/Admixture_mapping/test/output.txt
#SBATCH --error=/private/home/rsmerigl/codes/Admixture_mapping/test/error.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=512G

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/admixture
python /private/home/rsmerigl/codes/Admixture_mapping/test/pre_processing_testing.py