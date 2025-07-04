#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <result_folder> <data_folder>"
    exit 1
fi

result_folder="$1"
data_folder="$2"

ancestry_list=("AFR" "EAS" "EUR")

# Variables for the command
vcf_file="${data_folder}/input.vcf.gz"
# keep directory
keep_dir="${result_folder}/keep_files_chrom"
tmp_folder="${result_folder}/tmp"

keep_files=$(ls $keep_dir)

# File to save submitted job IDs

# loop on the keep files    
for keep_filename in $keep_files
do
    chrom=$(echo "$keep_filename" | awk -F'chr' '{print $2}' | awk -F'.' '{print $1}')
    snps_file="${result_folder}/snps_files/${keep_filename/keep/snps}"
    output_folder="${result_folder}/fine_mapping_ancestries_PCA_verbose/${keep_filename/_keep_chr$chrom.txt}_chr$chrom"
    mkdir -p $output_folder
    
    pheno=${keep_filename:4}
    pheno=${pheno/_keep_chr$chrom.txt/}
    ancestry=${keep_filename:0:3}

    phe_file="${data_folder}/phe_files/$pheno.phe"

    for ancestry_keep in ${ancestry_list[@]}
    do

        covar_file="${result_folder}/PCA_files/PCA_covar_files/${keep_filename/keep/covar}"
        covar_file="${covar_file%.txt}_${ancestry_keep}.tsv"

        output_file=$output_folder/${keep_filename/keep_chr$chrom.txt/output}
        output_file="${output_file:0:-6}chr${chrom}_output.${ancestry_keep}"

        keep_file="${result_folder}/keep_files_processed/${ancestry_keep}_keep.txt"

        if [ ! -f "$keep_file" ] || [ ! -f "$covar_file" ]; then
            echo "Keep file $keep_file or Covar file $covar_file does not exist. Skipping..."
            continue
        fi

        ../plink2 --vcf $vcf_file --pheno $phe_file --glm firth-fallback intercept --ci 0.95 --adjust --covar $covar_file --extract $snps_file --covar-variance-standardize --keep $keep_file --out $output_file --covar-col-nums 2-15

        glm_file="${output_file}.${pheno}.glm.logistic.hybrid"
        echo "Processing $glm_file for ancestry $ancestry_keep"

        if [ -f "$glm_file" ]; then
            echo "Filtering $glm_file (remove NA rows and assign ID)"

            tmp_filtered="${tmp_folder}/$(basename "$glm_file").filtered"

            awk 'NR == 1 || !/NA/' "$glm_file" > "$tmp_filtered"


            mv "$tmp_filtered" "$glm_file"
        fi




        done
    
done
