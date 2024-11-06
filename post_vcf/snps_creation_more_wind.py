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
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/{dataset}/snps_files_more_wind'
    wind_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/{dataset}/wind'
    vcf_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/{dataset}/vcf_file/ukbb.vcf.gz'
    tmp_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp'

    list_of_files = os.listdir(wind_folder)

    for wind_filename in list_of_files:  
        ancestry = wind_filename.split('_')[0]
        input_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/{dataset}/ancestry_{ancestry}.vcf'
        wind_file = os.path.join(wind_folder, wind_filename)
        wind_df_1 = pd.read_csv(wind_file, sep='\t')
        chroms = wind_df_1['chr'].unique()
        for chr in chroms:

            # Path for the output file
            output_file_1 = os.path.join(tmp_folder, 'chrom_pos_tmp.tsv')

            # Awk command to skip lines starting with '##' and extract only the CHROM and POS columns
            print(f'Starting awk command for chromosome {chr}')
            awk_command = f"awk '!/^##/ && $1 == \"{chr}\" {{print $1, $2}}' {input_file} > {output_file_1}"

            # Run the command using subprocess
            subprocess.run(awk_command, shell=True, check=True)

            # Load the output into a DataFrame
            df = pd.read_csv(output_file_1, sep=' ', header=None, names=['chr', 'pos'])

            df[['start', 'end']] = df['pos'].str.split('_', expand=True)

            # Convert 'start' and 'end' columns to integers
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)

            df = df.drop(columns=['pos'])

            # read the current window file into pandas and extract chrom, start, end
            wind_df = wind_df_1.loc[wind_df_1['chr'] == chr]

            n_wind = 3

            # Assume df and wind_df are sorted by 'chr' and 'start' (ascending)
            df = df.sort_values(['start']).reset_index(drop=True)
            wind_df = wind_df.sort_values(['start']).reset_index(drop=True)

            # Get the minimum start and maximum end from wind_df
            min_start = wind_df['start'].min()
            max_end = wind_df['end'].max()

            # Find rows in df that fall within this range
            extended_df = df[(df['start'] >= min_start) & (df['end'] <= max_end)].copy()

            # Get the index of the first and last rows in df that match min_start and max_end
            start_idx = df[df['start'] == min_start].index[0]
            end_idx = df[df['end'] == max_end].index[0]

            # Add the n_wind rows before the start and after the end
            start_extension = df.iloc[max(0, start_idx - n_wind):start_idx]
            end_extension = df.iloc[end_idx + 1:end_idx + 1 + n_wind]

            # Concatenate the extensions with the main range
            extended_wind_df = pd.concat([start_extension, extended_df, end_extension]).drop_duplicates().reset_index(drop=True)
            
            output_file = os.path.join(output_folder, wind_filename.replace('wind', f'snps_chr{chr}'))
            tmp_file_vcf = os.path.join(tmp_folder, wind_filename.replace('wind.txt', 'tmp.vcf'))
            tmp_file_snps = os.path.join(tmp_folder, wind_filename.replace('wind', 'tmp_snps'))

            # define the snps list
            list_of_snps = []

            # read the current window file into pandas and extract chrom, start, end
            wind_df = wind_df_1.loc[wind_df_1['chr'] == chr]
            for index, row in extended_wind_df.iterrows():
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

            last_df = pd.DataFrame(list_of_snps, columns=['#ID'])
            last_df.to_csv(output_file, index=False)
            print(f"File {output_file} chr {chrom} created")




        