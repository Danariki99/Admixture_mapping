import sys
import subprocess
import os
import shutil
import pandas as pd


if __name__ == '__main__':
    # Check if the dataset argument is provided
    if len(sys.argv) != 2:
        print("Usage: python post_processing.py <dataset>")
        sys.exit(1)

    # The dataset variable is taken from the command line argument
    dataset = sys.argv[1]

    # Use the dataset variable to construct file paths
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/{dataset}/snps_files'
    wind_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/{dataset}/wind'
    vcf_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/{dataset}/vcf_file/ukbb.vcf.gz'
    tmp_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp'

    list_of_files = os.listdir(wind_folder)

    for wind_filename in list_of_files:  # Process just one file for demonstration
        wind_file = os.path.join(wind_folder, wind_filename)
        output_file = os.path.join(output_folder, wind_filename.replace('wind', 'snps'))
        tmp_file_vcf = os.path.join(tmp_folder, wind_filename.replace('wind.txt', 'tmp.vcf'))
        tmp_file_snps = os.path.join(tmp_folder, wind_filename.replace('wind', 'tmp_snps'))

        # define the snps list
        list_of_snps = []

        # read the current window file into pandas and extract chrom, start, end
        wind_df = pd.read_csv(wind_file, sep='\t')
        for index, row in wind_df.iterrows():
            chrom = row['chr']
            start = row['start']
            end = row['end']

            # Create a temporary vcf file with only the SNPs in the window
            plink_command = [
                "/private/home/rsmerigl/plink2",    
                "--vcf", vcf_file,                  
                "--chr", str(chrom),               
                "--from-bp", str(start),            
                "--to-bp", str(end),                
                "--recode", "vcf",                  
                "--out", tmp_file_vcf.replace('.vcf','')                   
            ]

            subprocess.run(plink_command)

            #extract only the variants ID in a tmp file
            awk_command = f"awk '!/^##/ {{print $3}}' {tmp_file_vcf} > {tmp_file_snps}"

            subprocess.run(awk_command, shell=True, check=True)

            #read the file in pandas and extract the list of snps
            df = pd.read_csv(tmp_file_snps)
            list_of_snps.extend(df['ID'].values)

            # delete the folder and recreate it to free space
            shutil.rmtree(tmp_folder)
            os.makedirs(tmp_folder)

        last_df = pd.DataFrame(list_of_snps, columns=['ID'])
        last_df.to_csv(output_file, index=False)
        print(f"File {output_file} created")




        