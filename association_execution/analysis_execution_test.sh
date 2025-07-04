#!/bin/bash


# Stop execution if any command fails
set -e

# Check if an argument is passed
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <results_folder> <data_dir>"
  exit 1
fi

# Input parameter (e.g. dataset name or base folder)
output_root="$1"
data_dir="$2"

# Define base directories
vcf_dir="$output_root/vcf_files"
pheno_dir="$data_dir/phe_files"
covar_file="$data_dir/input.covar"
keep_file="$data_dir/input.keep"

# Loop over each VCF file
for vcf_file in "$vcf_dir"/*.vcf; do
    # Extract ancestry from filename: e.g., ancestry_EUR.vcf -> EUR
    filename=$(basename "$vcf_file")
    ancestry=$(echo "$filename" | sed -E 's/^ancestry_([A-Z]+)\.vcf$/\1/')

    # Loop over each phenotype file
    for pheno_file in "$pheno_dir"/*.phe; do
        pheno_base=$(basename "$pheno_file" .phe)
        echo "Running PLINK2 for ancestry $ancestry and phenotype $pheno_base"

        # Define output directory and file base name
        output_dir="$output_root/output_ancestry_${ancestry}/${pheno_base}"
        mkdir -p "$output_dir"
        output_file="$output_dir/output"

        # Run plink2
        ../plink2 \
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
