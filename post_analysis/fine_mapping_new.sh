#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset>"
    exit 1
fi

dataset=$1

vcf_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/vcf_file/ukbb.vcf.gz"
wind_folder="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/${dataset}/wind"
covar_folder="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new"
keep_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id"
base_output="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/fine_mapping_new"
sbatch_dir="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/${dataset}/fine_mapping_new"

mkdir -p $sbatch_dir

job_ids_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp/submitted_job_ids_fine_mapping_new.txt"
> $job_ids_file

job_counter=1

for wind_file in ${wind_folder}/*.txt; do
    wind_filename=$(basename $wind_file)
    hit_ancestry=$(echo $wind_filename | cut -d'_' -f1)
    pheno=$(echo $wind_filename | cut -d'_' -f2)

    phe_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/${dataset}/${pheno}.phe"

    # Read each significant window (skip header)
    while IFS=$'\t' read -r chr start end rest; do

        output_folder="${base_output}/${hit_ancestry}_${pheno}_chr${chr}"
        mkdir -p $output_folder

        # Per-window covariate file: base covariates + LAI of this specific window
        # Columns: IID(1) age(2) sex(3) BMI(4) AFR(5) AHG(6) EAS(7) EUR(8) NAT(9) OCE(10) SAS(11) WAS(12) LAI(13)
        # EUR skipped (col 8) as reference category, same as original admixture mapping
        covar_file="${covar_folder}/${hit_ancestry}_${pheno}_chr${chr}_${start}_${end}_covar.tsv"

        output_file="${output_folder}/${hit_ancestry}_${pheno}_chr${chr}_${start}_${end}"

        command_to_run="/private/home/rsmerigl/plink2 --vcf $vcf_file --pheno $phe_file --glm firth-fallback --ci 0.95 --adjust --covar $covar_file --chr $chr --from-bp $start --to-bp $end --covar-variance-standardize --keep $keep_file --out $output_file --covar-col-nums 2-7,9-13"

        sbatch_file="${sbatch_dir}/job_${hit_ancestry}_${pheno}_chr${chr}_${start}_${end}_${job_counter}.sh"

        echo "#!/bin/bash" > $sbatch_file
        echo "#SBATCH --partition=long" >> $sbatch_file
        echo "#SBATCH --job-name=fm_${hit_ancestry}_${pheno}_${chr}" >> $sbatch_file
        echo "#SBATCH --output=${sbatch_file%.sh}.out" >> $sbatch_file
        echo "#SBATCH --error=${sbatch_file%.sh}.err" >> $sbatch_file
        echo "#SBATCH --nodes=1" >> $sbatch_file
        echo "#SBATCH --cpus-per-task=1" >> $sbatch_file
        echo "#SBATCH --time=1-00:00:00" >> $sbatch_file
        echo "#SBATCH --mem=128G" >> $sbatch_file
        echo "$command_to_run" >> $sbatch_file

        sbatch_output=$(sbatch $sbatch_file)
        job_id=$(echo $sbatch_output | grep -oP '(?<=Submitted batch job )\d+')

        echo "Submitted job ID: $job_id ($hit_ancestry $pheno chr$chr $start-$end)"
        echo "$job_id" >> $job_ids_file

        ((job_counter++))

    done < <(tail -n +2 $wind_file)

done
