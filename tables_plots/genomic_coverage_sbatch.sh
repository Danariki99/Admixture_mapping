#!/bin/bash
#SBATCH --job-name=genomic_coverage
#SBATCH --error=./genomic_coverage.err
#SBATCH --output=./genomic_coverage.out
#SBATCH --partition=long
#SBATCH --mem=64G
#SBATCH --time=12:00:00

source /private/home/rsmerigl/anaconda3/etc/profile.d/conda.sh
conda activate r_env

python /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/genomic_coverage_plot.py
