#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=single_ancestry
#SBATCH --output=/private/home/rsmerigl/codes/Admixture_mapping/test/output_merging_single_ancestry.txt
#SBATCH --error=/private/home/rsmerigl/codes/Admixture_mapping/test/error_merging_single_ancestry.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/admixture
python /private/home/rsmerigl/codes/Admixture_mapping/test/merging_single_ancestry.py