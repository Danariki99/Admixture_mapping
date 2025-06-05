#!/bin/bash

# Stop execution if any command fails
set -e

# Define base directories
vcf_dir="./data/vcf_files"
pheno_dir="./data/phe_files"
covar_file="./data/input.covar"
keep_file="./data/input.keep"
output_root="./results"

# Loop over each VCF file
for vcf_file in "$vcf_dir"/*.vcf; do
    # Extract ancestry from filename: e.g., chr1_AFR.vcf -> AFR
    filename=$(basename "$vcf_file")
    ancestry=$(echo "$filename" | sed -E 's/^.*_([A-Z]+)\.vcf$/\1/')

    # Loop over each phenotype file
    for pheno_file in "$pheno_dir"/*.phe; do
        pheno_base=$(basename "$pheno_file" .phe)
        echo "Running PLINK2 for ancestry $ancestry and phenotype $pheno_base"

        # Define output directory and file base name
        output_dir="$output_root/output_ancestry_${ancestry}/${pheno_base}"
        mkdir -p "$output_dir"
        output_file="$output_dir/output"

        # Run plink2
        /private/home/rsmerigl/plink2 \
            --vcf "$vcf_file" \
            --pheno "$pheno_file" \
            --glm firth-fallback hide-covar \
            --ci 0.95 \
            --adjust \
            --covar "$covar_file" \
            --covar-variance-standardize \
            --keep "$keep_file" \
            --out "$output_file" \
            --covar-col-nums 2-14

    done
done
