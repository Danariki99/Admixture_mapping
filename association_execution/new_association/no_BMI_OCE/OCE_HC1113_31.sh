#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=job_OCE_HC1113_31
#SBATCH --output=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/new_association/no_BMI_OCE/OCE_HC1113_31.out
#SBATCH --error=/private/home/rsmerigl/codes/Admixture_mapping/association_execution/new_association/no_BMI_OCE/OCE_HC1113_31.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G
/private/home/rsmerigl/plink2 --vcf /private/groups/ioannidislab/smeriglio/merged_vcfs/ancestry_OCE.vcf --pheno /private/groups/ioannidislab/smeriglio/phe_files/HC1113.phe --glm firth-fallback hide-covar --ci 0.95 --adjust --covar /private/groups/ioannidislab/smeriglio/covar_file/ukb24983_GWAS_covar_new.phe --covar-variance-standardize --keep /private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt --out /private/groups/ioannidislab/smeriglio/output_OCE/no_BMI/HC1113/output --covar-col-nums 2,3,93
