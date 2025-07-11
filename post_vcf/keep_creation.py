import sys
import subprocess
import os
import pandas as pd

# This script is used to create keep files with only the samples that, in the selected windows,
# are fully of the detected ancestry. The script reads the windows files, filters the VCF files
# to include only the selected windows, and then filters the samples to include only the samples
# that are fully of the detected ancestry in the selected windows. The script also filters the
# samples to include only the samples that are in the old keep file.

if __name__ == '__main__':
    # Check if the dataset argument is provided
    if len(sys.argv) != 2:
        print("Usage: python post_processing.py <dataset>")
        sys.exit(1)

    # The dataset variable is taken from the command line argument
    dataset = sys.argv[1]

    # Use the dataset variable to construct file paths
    keep_output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/{dataset}/keep_files'
    wind_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/{dataset}/wind'
    vcf_generic_filename = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/{dataset}/ancestry_*.vcf'
    tmp_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp'
    old_keep_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/{dataset}/keep_file.txt'

    list_of_files = os.listdir(wind_folder)

    for wind_filename in list_of_files:  # Process just one file for demonstration
        wind_file = os.path.join(wind_folder, wind_filename)
        output_file = os.path.join(keep_output_folder, wind_filename.replace('wind', 'keep'))
        tmp_file = os.path.join(tmp_folder, wind_filename.replace('wind', 'tmp'))

        # Assuming the VCF file is named similarly to the wind file, construct the VCF file path
        vcf_file = vcf_generic_filename.replace('*', wind_filename.split('_')[0])

        # read the old keep file
        old_keep = pd.read_csv(old_keep_file)

        # Construct the awk command to filter windows and include #CHROM
        print(f'Filtering windows file {wind_filename}...')
        awk_command = f"""
        awk 'NR==FNR {{win[$1 FS $2 "_" $3]; next}} 
        FNR==1 {{if ($0 ~ /^#CHROM/) print $0; next}}  # Print only #CHROM line as header
        $1 ~ /^#CHROM/ {{print; next}}  # Print the line that starts with #CHROM
        ($1 FS $2) in win {{print}}' {wind_file} {vcf_file}
        """

        # Execute the awk command to generate the temporary VCF
        with open(tmp_file, 'w') as tmp_f:
            subprocess.run(awk_command, shell=True, stdout=tmp_f)

        # Now use awk to filter samples with only 1/1 values directly from tmp_file
        print('Filtering samples...')
        awk_filter_command = f"""
        awk 'NR == 1 {{ # Process the first line (header)
            for (i = 10; i <= NF; i++) {{
                sample[i] = $i;  # Store sample names
                all_one_one[i] = 1;  # Assume all are 1/1 initially
            }}
            next;
        }} 
        {{  
            for (i = 10; i <= NF; i++) {{
                if ($i != "1/1") {{
                    all_one_one[i] = 0;  # Mark as not all 1/1 if any row is different
                }}
            }}
        }} 
        END {{
            print "#IID";
            for (i = 10; i <= NF; i++) {{
                if (all_one_one[i]) {{
                    print sample[i];  # Print sample only if all rows are 1/1
                }}
            }}
        }}' {tmp_file}
    """

        # Write the selected samples to the output file
        with open(output_file, 'w') as out_f:
            subprocess.run(awk_filter_command, shell=True, stdout=out_f)

        # Read the selected samples from the output file
        selected_samples = pd.read_csv(output_file)

        # Filter the new keep file to include only samples that are in the old keep file
        new_keep = selected_samples[selected_samples['#IID'].isin(old_keep['#IID'])]

        # Write the new keep file
        new_keep.to_csv(output_file, index=False)

        print(f'Finished processing {wind_filename}\nThe number of samples in the new keep file is {len(new_keep)}')

        # Clean up the temporary files
        os.remove(tmp_file)




        
