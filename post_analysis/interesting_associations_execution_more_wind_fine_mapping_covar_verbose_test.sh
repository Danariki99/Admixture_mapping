#!/bin/bash

# Stop if any command fails
set -e

# === Controllo input ===
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <vcf_file.vcf.gz>"
    exit 1
fi

vcf_file="$1"

# === Configurazione pipeline ===

# Define ancestry list
ancestry_list=("AFR" "AHG" "EAS" "EUR" "NAT" "SAS" "WAS")

# Percorsi
phe_folder="./data/phe_files"
keep_folder="./results/keep_files_processed"
snps_folder="./results/snps_files"
covar_folder="./results/pca/PCA_covar_files"
output_root="./results/fine_mapping_verbose"
tmp_folder="./results/tmp"

mkdir -p "$output_root"
mkdir -p "$tmp_folder"

# Loop su tutti i file .snps
for snps_file in "$snps_folder"/*.snps; do
    filename=$(basename "$snps_file")
    echo "Processing $filename"

    # Estrae PHENO e CHR dal nome file: <PHENO>_<ANCESTRY>_chrN.snps
    pheno=$(echo "$filename" | cut -d'_' -f2)
    chr=$(echo "$filename" | grep -oP 'chr\d+')

    phe_file="$phe_folder/${pheno}.phe"

    for ancestry_keep in "${ancestry_list[@]}"; do
        keep_file="$keep_folder/${ancestry_keep}_keep.txt"
        covar_file="${covar_folder}/covar_${pheno}_${ancestry_keep}.tsv"

        output_dir="${output_root}/${pheno}/${chr}/${ancestry_keep}"
        mkdir -p "$output_dir"
        output_file="${output_dir}/output"

        echo "→ Running PLINK2 for $pheno, $ancestry_keep, $chr"

        /private/home/rsmerigl/plink2 \
            --vcf "$vcf_file" \
            --pheno "$phe_file" \
            --glm firth-fallback intercept \
            --ci 0.95 \
            --adjust \
            --covar "$covar_file" \
            --extract "$snps_file" \
            --covar-variance-standardize \
            --keep "$keep_file" \
            --out "$output_file" \
            --covar-col-nums 2-15
    done
done
