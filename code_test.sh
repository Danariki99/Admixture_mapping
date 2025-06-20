#!/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <data_folder> <result_folder>"
    exit 1
fi

# Assign the arguments to variables
data_folder="$1"
result_folder="$2"

vcf_folder=$(python LAI/chrom_division.py --results "$result_folder" --data "$data_folder")

echo "Using GnoMix software"
python LAI/gnomix_training_test.py --vcf_folder "$vcf_folder" --panel_folder "$panel_folder"

msp_folder=$(./LAI/files_moving.sh "$result_folder")


# Pass the dataset variable to the Python scripts
python pre_processing/pre_processing_testy.py "$msp_folder" "$result_folder"

# Execute job_submission_BMI.sh and capture the path of the file containing job IDs
./association_execution/analysis_execution_test.sh  "$result_folder" "$data_folder"

python post_processing/post_processing_test.py "$result_folder" "$data_folder"

python post_analysis/ancestry_counts_test.py "$msp_folder" "$result_folder"

python post_analysis/keep_files_test.py "$result_folder"    

python post_analysis/keep_files_post_processing_test.py "$result_folder" "$data_folder"

python post_vcf/snps_creation_more_wind_test.py "$result_folder" "$data_folder"

python post_analysis/LD_pruning_test.py "$result_folder" "$data_folder"

python post_analysis/PCA_generation_test.py "$result_folder" "$data_folder"

python post_analysis/new_covariable_creation_test.py "$result_folder" "$data_folder" 

python post_analysis/PCA_covar_files_creation_test.py "$result_folder"

python post_vcf/keep_creation_chrom_test.py "$result_folder" "$data_folder"

./post_analysis/interesting_associations_execution_more_wind_fine_mapping_covar_verbose_test.sh "$result_folder" "$data_folder"

python post_analysis/betas_extraction_test.py "$result_folder"

python post_analysis/probabilities_samples_dataset_extraction_test.py "$result_folder"

python post_analysis/delta_probabilities_computation_and_plot_v2_test.py "$result_folder"





