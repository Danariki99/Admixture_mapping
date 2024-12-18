#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset>"
    exit 1
fi

# Assign the first argument to the variable 'dataset'
dataset=$1

# Directory where sbatch scripts will be created
sbatch_dir="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/${dataset}"

# Base name for sbatch scripts
sbatch_base="job"

# Create the sbatch scripts directory if it doesn't exist
mkdir -p $sbatch_dir

# Variables for the command
vcf_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/vcf_file/ukbb.vcf.gz"
# keep directory
keep_dir="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/keep_files_chrom"

keep_files=$(ls $keep_dir)

# File to save submitted job IDs
job_ids_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp/submitted_job_ids.txt"
# Clear the file to start fresh
> $job_ids_file

# Create and submit sbatch scripts
job_counter=1

# loop on the keep files    
for keep_filename in $keep_files
do
    chrom=$(echo "$keep_filename" | awk -F'chr' '{print $2}' | awk -F'.' '{print $1}')
    snps_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/snps_files_more_wind_3/${keep_filename/keep/snps}"
    output_folder=/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/${dataset}/fine_mapping_ancestries/${keep_filename/_keep_chr$chrom.txt}_chr$chrom
    mkdir -p $output_folder

    output_file=$output_folder/${keep_filename/keep_chr$chrom.txt/output}
    output_file="${output_file:0:-6}chr${chrom}_output"

    
    pheno=${keep_filename:4}
    pheno=${pheno/_keep_chr$chrom.txt/}
    ancestry=${keep_filename:0:3}

    covar_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files/${keep_filename/keep/covar}"

    keep_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/${dataset}/keep_files/${ancestry}_keep.txt"

    phe_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/${dataset}/$pheno.phe"

    command_to_run="/private/home/rsmerigl/plink2 --vcf $vcf_file --pheno $phe_file --glm firth-fallback hide-covar --ci 0.95 --adjust --covar $covar_file --extract $snps_file --covar-variance-standardize --keep $keep_file --out $output_file --covar-col-nums 2-15"

    sbatch_file="$sbatch_dir/${keep_filename/keep.txt/snps}_$job_counter.sh"

    echo "#!/bin/bash" > $sbatch_file
    echo "#SBATCH --partition=long" >> $sbatch_file
    echo "#SBATCH --job-name=${sbatch_base}_${keep_filename/keep.txt/snps}_$job_counter" >> $sbatch_file
    echo "#SBATCH --output=$sbatch_dir/${keep_filename/keep.txt/snps}_$job_counter.out" >> $sbatch_file
    echo "#SBATCH --error=$sbatch_dir/${keep_filename/keep.txt/snps}_$job_counter.err" >> $sbatch_file
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
