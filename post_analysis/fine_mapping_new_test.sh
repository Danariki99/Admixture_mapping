#!/bin/bash
# TEST version of fine_mapping_new.sh — CONDITIONAL fine mapping.
# Same logic as production but runs locally on the data we pass (no SLURM), using
# result_folder / data_folder. For each significant window it runs plink2 --glm
# with the per-window covariate file that includes the LAI (col 13), so the model
# yields both the SNP (ADD) and the local-ancestry (LAI) test rows.

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <result_folder> <data_folder>"
    exit 1
fi

result_folder="$1"
data_folder="$2"

vcf_file="${data_folder}/input.vcf.gz"
wind_folder="${result_folder}/FUMA/wind"
covar_folder="${result_folder}/wind_covar_files_new"
base_output="${result_folder}/fine_mapping_new"
# optional keep file (all unrelated samples); if absent, all samples are used
keep_file="${result_folder}/kcutoff/kcutoff.king.cutoff.in.id"

mkdir -p "$base_output"

keep_arg=""
[ -f "$keep_file" ] && keep_arg="--keep $keep_file"

for wind_file in ${wind_folder}/*.txt; do
    [ -e "$wind_file" ] || continue
    wind_filename=$(basename "$wind_file")
    hit_ancestry=$(echo "$wind_filename" | cut -d'_' -f1)
    pheno=$(echo "$wind_filename" | cut -d'_' -f2)

    phe_file="${data_folder}/phe_files/${pheno}.phe"

    # read each significant window (skip header): chr start end ...
    while IFS=$'\t' read -r chr start end rest; do

        output_folder="${base_output}/${hit_ancestry}_${pheno}_chr${chr}"
        mkdir -p "$output_folder"

        # per-window covariate file: base covariates + LAI of THIS window
        # cols: IID(1) age(2) sex(3) BMI(4) AFR(5) AHG(6) EAS(7) EUR(8) NAT(9) OCE(10) SAS(11) WAS(12) LAI(13)
        # EUR (col 8) dropped as reference, LAI (col 13) included -> ADD + LAI tests
        covar_file="${covar_folder}/${hit_ancestry}_${pheno}_chr${chr}_${start}_${end}_covar.tsv"
        output_file="${output_folder}/${hit_ancestry}_${pheno}_chr${chr}_${start}_${end}"

        if [ ! -f "$covar_file" ]; then
            echo "Covar file $covar_file missing, skipping window $chr:$start-$end"
            continue
        fi

        ../plink2 --vcf "$vcf_file" --pheno "$phe_file" --glm firth-fallback --ci 0.95 \
            --adjust --covar "$covar_file" --chr "$chr" --from-bp "$start" --to-bp "$end" \
            --covar-variance-standardize $keep_arg --out "$output_file" --covar-col-nums 2-7,9-13

    done < <(tail -n +2 "$wind_file")

done
