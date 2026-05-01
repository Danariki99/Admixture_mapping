#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset>"
    exit 1
fi

dataset=$1

vcf_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/vcf_file/ukbb.vcf.gz"
covar_folder="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new_6wind"
keep_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id"
base_output="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/fine_mapping_new_6wind"
sbatch_dir="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/${dataset}/fine_mapping_new_6wind"

mkdir -p $sbatch_dir

job_ids_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp/submitted_job_ids_fine_mapping_new_6wind.txt"
> $job_ids_file

job_counter=1

sig_covar_folder="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new"

# Iterate only over extension covariate files (sig windows already processed in fine_mapping_new)
for covar_file in ${covar_folder}/*.tsv; do
    # Skip if this window was already processed in the sig-only run
    if [ -f "${sig_covar_folder}/$(basename $covar_file)" ]; then
        continue
    fi
    covar_filename=$(basename $covar_file)
    # Filename: {ancestry}_{pheno}_chr{N}_{start}_{end}_covar.tsv
    base=${covar_filename/_covar.tsv/}   # AFR_HC1574_chr5_87070639_87391401
    parts=(${base//_/ })
    hit_ancestry=${parts[0]}
    pheno=${parts[1]}
    chr_num=${parts[2]//chr/}
    start=${parts[3]}
    end=${parts[4]}

    phe_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/${dataset}/${pheno}.phe"

    output_folder="${base_output}/${hit_ancestry}_${pheno}_chr${chr_num}"
    mkdir -p $output_folder

    output_file="${output_folder}/${hit_ancestry}_${pheno}_chr${chr_num}_${start}_${end}"

    command_to_run="/private/home/rsmerigl/plink2 --vcf $vcf_file --pheno $phe_file --glm firth-fallback --ci 0.95 --adjust --covar $covar_file --chr $chr_num --from-bp $start --to-bp $end --covar-variance-standardize --keep $keep_file --out $output_file --covar-col-nums 2-7,9-13"

    sbatch_file="${sbatch_dir}/job_${hit_ancestry}_${pheno}_chr${chr_num}_${start}_${end}_${job_counter}.sh"

    echo "#!/bin/bash" > $sbatch_file
    echo "#SBATCH --partition=long" >> $sbatch_file
    echo "#SBATCH --job-name=fm6_${hit_ancestry}_${pheno}_${chr_num}" >> $sbatch_file
    echo "#SBATCH --output=${sbatch_file%.sh}.out" >> $sbatch_file
    echo "#SBATCH --error=${sbatch_file%.sh}.err" >> $sbatch_file
    echo "#SBATCH --nodes=1" >> $sbatch_file
    echo "#SBATCH --cpus-per-task=1" >> $sbatch_file
    echo "#SBATCH --time=1-00:00:00" >> $sbatch_file
    echo "#SBATCH --mem=128G" >> $sbatch_file
    echo "$command_to_run" >> $sbatch_file

    sbatch_output=$(sbatch $sbatch_file)
    job_id=$(echo $sbatch_output | grep -oP '(?<=Submitted batch job )\d+')

    echo "Submitted job ID: $job_id ($hit_ancestry $pheno chr$chr_num $start-$end)"
    echo "$job_id" >> $job_ids_file

    ((job_counter++))

done
