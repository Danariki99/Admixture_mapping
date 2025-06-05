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
