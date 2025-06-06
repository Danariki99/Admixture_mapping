#!/bin/bash

# Check if a dataset argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <msp_folder>"
    exit 1
fi

# Assign the first argument to the variable 'dataset'
file=$1

# Pass the dataset variable to the Python scripts
python pre_processing/pre_processing_testy.py "$file"

# Execute job_submission_BMI.sh and capture the path of the file containing job IDs
./association_execution/analysis_execution_test.sh 

python post_processing/post_processing_test.py

python post_analysis/ancestry_counts_test.py "$file"

python post_analysis/keep_files_test.py

python post_analysis/keep_files_post_processing_test.py

python post_vcf/snps_creation_more_wind_test.py

python post_analysis/PCA_generation_test.py

./post_analysis/interesting_associations_execution_more_wind_fine_mapping_covar_verbose_test.sh

python post_analysis/betas_extraction_test.py

python post_analysis/probabilities_samples_dataset_extraction_test.py

python post_analysis/delta_probabilities_computation_and_plot_v2_test.py





