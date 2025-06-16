#!/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <vcf_file> <panel_folder>"
    exit 1
fi

# Assign the arguments to variables
file=$1
panel_folder=$2


vcf_folder=$(python LAI/chrom_division.py --vcf "$file")

echo "Using GnoMix software"
python LAI/gnomix_training_test.py --vcf_folder "$vcf_folder" --panel_folder "$panel_folder"

msp_folder=$(python LAI/gnomix.py --vcf_folder "$vcf_folder" --gnomix_ref "$panel_folder")

msp_folder=$(./LAI/files_moving.sh)


# Pass the dataset variable to the Python scripts
python pre_processing/pre_processing_testy.py "$msp_folder"

# Execute job_submission_BMI.sh and capture the path of the file containing job IDs
./association_execution/analysis_execution_test.sh 

python post_processing/post_processing_test.py

python post_analysis/ancestry_counts_test.py "$msp_folder"

python post_analysis/keep_files_test.py

python post_analysis/keep_files_post_processing_test.py

python post_vcf/snps_creation_more_wind_test.py

python post_analysis/LD_pruning_test.py

python post_analysis/PCA_generation_test.py "$file"

./post_analysis/interesting_associations_execution_more_wind_fine_mapping_covar_verbose_test.sh "$file"

python post_analysis/betas_extraction_test.py

python post_analysis/probabilities_samples_dataset_extraction_test.py

python post_analysis/delta_probabilities_computation_and_plot_v2_test.py





