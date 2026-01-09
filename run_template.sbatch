#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=serics_gpu
#SBATCH --job-name=admixture
#SBATCH --output=../slurm_outputs/admixture.out
#SBATCH --error=../slurm_outputs/admixture.err
#SBATCH --mem=500G

singularity run --bind ../data_results_admixture_test/:/src singularity.sif /src/data /src/results
