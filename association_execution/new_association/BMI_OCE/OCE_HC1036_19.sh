#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=job_OCE_HC1036_19
#SBATCH --output=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/new_association/BMI_OCE/OCE_HC1036_19.out
#SBATCH --error=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/new_association/BMI_OCE/OCE_HC1036_19.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G
/private/home/rsmerigl/plink2 --vcf /private/groups/ioannidislab/smeriglio/merged_vcfs/ancestry_OCE.vcf --pheno /private/groups/ioannidislab/smeriglio/phe_files/HC1036.phe --glm firth-fallback hide-covar --ci 0.95 --adjust --covar /private/groups/ioannidislab/smeriglio/covar_file/ukb24983_GWAS_covar_new.phe --covar-variance-standardize --keep /private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt --out /private/groups/ioannidislab/smeriglio/output_OCE/BMI/HC1036/output --covar-col-nums 2-4,93
