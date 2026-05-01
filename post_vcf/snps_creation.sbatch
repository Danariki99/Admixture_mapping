#!/bin/bash
#SBATCH --job-name=snps_creation_more_wind
#SBATCH --output=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_vcf/snps_creation.out
#SBATCH --error=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_vcf/snps_creation.err
#SBATCH --partition=long
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

source /private/home/rsmerigl/anaconda3/bin/activate test

cd /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_vcf

python -u snps_creation_more_wind.py ukbb
