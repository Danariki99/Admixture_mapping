#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=pipeline_execution
#SBATCH --output=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/execution_output.txt
#SBATCH --error=/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/execution_error.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=512G

# Check if a dataset argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset>"
    exit 1
fi

# Assign the first argument to the variable 'dataset'
dataset=$1

source /private/home/rsmerigl/anaconda3/bin/activate /private/home/rsmerigl/anaconda3/envs/r_env

# Pass the dataset variable to the Python scripts
python pre_processing/pre_processing.py $dataset

# Execute job_submission_BMI.sh and capture the path of the file containing job IDs
./association_execution/job_submission_BMI.sh $dataset

job_ids_file="/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp/submitted_job_ids.txt"

# Read job IDs from the file
job_ids=$(cat "$job_ids_file")

# Loop until all jobs are completed
all_done=0
while [ $all_done -eq 0 ]; do
    all_done=1
    for job_id in $job_ids; do
        # Check if the job is still in the queue or running
        if squeue --job $job_id 2>&1 | grep -q "$job_id"; then
            all_done=0
            break
        fi
    done
    
    # If not all jobs are completed, wait an hour before checking again
    if [ $all_done -eq 0 ]; then
        sleep 1h
    fi
done
echo "All jobs are completed, starting post processing"
# All jobs are completed, execute the desired command
python post_processing/post_processing.py $dataset

# Cleanup: remove the job IDs file
rm "$job_ids_file"