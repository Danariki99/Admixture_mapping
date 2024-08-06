#!/bin/bash

# Directory where sbatch scripts will be created
sbatch_dir="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/ukbb"

# Base name for sbatch scripts
sbatch_base="job"

# Ancestry array
ancestry=("AFR" "AHG" "EAS" "EUR" "NAT" "SAS" "WAS")

# Get list of files in pheno folder
pheno_folder="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/ukbb"
pheno_files=($(ls $pheno_folder))

# Create the sbatch scripts directory if it doesn't exist
mkdir -p $sbatch_dir

# Variables for the command
covar_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered.phe"
keep_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/ukbb/keep_file.txt"

# File to save submitted job IDs
job_ids_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp/submitted_job_ids.txt"
# Clear the file to start fresh
> $job_ids_file

# Create and submit sbatch scripts
job_counter=1
for anc in "${ancestry[@]}"
do
    vcf_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/ukbb/ancestry_${anc}.vcf"

    for pheno_file in "${pheno_files[@]}"
    do
        current_pheno_file="$pheno_folder/$pheno_file"
        
        pheno_base=$(basename $pheno_file)
        pheno_base=${pheno_base%.phe}
        echo "Processing $pheno_base"
                
        output_dir="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb/output_ancestry_${anc}/$pheno_base"
        mkdir -p $output_dir
        output_file="$output_dir/output"

        command_to_run="/private/home/rsmerigl/plink2 --vcf $vcf_file --pheno $current_pheno_file --glm firth-fallback hide-covar --ci 0.95 --adjust --covar $covar_file --covar-variance-standardize --keep $keep_file --out $output_file --covar-col-nums 2-14"

        sbatch_file="$sbatch_dir/${anc}_${pheno_base}_$job_counter.sh"
        echo "#!/bin/bash" > $sbatch_file
        echo "#SBATCH --partition=long" >> $sbatch_file
        echo "#SBATCH --job-name=${sbatch_base}_${anc}_${pheno_base}_$job_counter" >> $sbatch_file
        echo "#SBATCH --output=$sbatch_dir/${anc}_${pheno_base}_$job_counter.out" >> $sbatch_file
        echo "#SBATCH --error=$sbatch_dir/${anc}_${pheno_base}_$job_counter.err" >> $sbatch_file
        echo "#SBATCH --nodes=1" >> $sbatch_file
        echo "#SBATCH --cpus-per-task=1" >> $sbatch_file
        echo "#SBATCH --time=7-00:00:00" >> $sbatch_file
        echo "#SBATCH --mem=512G" >> $sbatch_file
        echo "$command_to_run" >> $sbatch_file

        # Submit sbatch script and capture the output
        sbatch_output=$(sbatch $sbatch_file)
        # Extract job ID
        job_id=$(echo $sbatch_output | grep -oP '(?<=Submitted batch job )\d+')
        echo "Submitted job ID: $job_id"

        # Save the job ID to the file
        echo "$job_id" >> $job_ids_file

        ((job_counter++))
    done
done