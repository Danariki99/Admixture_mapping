#!/bin/bash

# Directory where sbatch scripts will be created
sbatch_dir="/private/home/rsmerigl/codes/Admixture_mapping/association_execution/new_association/no_BMI"

# Base name for sbatch scripts
sbatch_base="job"

# Ancestry array
#ancestry=("AFR" "AHG" "EAS" "EUR" "NAT" "OCE" "SAS" "WAS")
ancestry=("OCE" "SAS" "WAS")

# Get list of files in pheno folder
pheno_folder="/private/groups/ioannidislab/smeriglio/phe_files"
pheno_files=($(ls $pheno_folder))

# Create the sbatch scripts directory if it doesn't exist
mkdir -p $sbatch_dir

# Variables for the command
covar_file="/private/groups/ioannidislab/smeriglio/covar_file/ukb24983_GWAS_covar.phe"
keep_file="/private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt"

# Create and submit sbatch scripts
job_counter=546
for anc in "${ancestry[@]}"
do
    vcf_file="/private/groups/ioannidislab/smeriglio/merged_vcfs/ancestry_${anc}.vcf"

    for pheno_file in "${pheno_files[@]}"
    do
        # Update the current_pheno_file and output_file variables for each combination of ancestry and phenotype
        current_pheno_file="$pheno_folder/$pheno_file"
        
        # Get the base name of the phenotype file
        pheno_base=$(basename $pheno_file)
        pheno_base=${pheno_base%.phe}
        echo "Processing $pheno_base"
                
        # Construct the output file path
        output_file="/private/groups/ioannidislab/smeriglio/output/no_BMI/output_ancestry_${anc}/$pheno_base/output"

        # Command to run
        command_to_run="/private/home/rsmerigl/plink2 --vcf $vcf_file --pheno $current_pheno_file --glm firth-fallback hide-covar --ci 0.95 --adjust --covar $covar_file --covar-variance-standardize --keep $keep_file --out $output_file --covar-col-nums 2,3,48-57"

        # Create sbatch script
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


        # Submit sbatch script
        sbatch $sbatch_file

        # Increment job counter
        ((job_counter++))
    done
done