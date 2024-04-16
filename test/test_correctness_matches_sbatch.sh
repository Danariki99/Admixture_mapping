#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=correctness_matches
#SBATCH --output=/private/home/rsmerigl/codes/Admixture_mapping/test/output_matches.txt
#SBATCH --error=/private/home/rsmerigl/codes/Admixture_mapping/test/error_matches.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/admixture
python /private/home/rsmerigl/codes/Admixture_mapping/test/test_correctness_matches.py