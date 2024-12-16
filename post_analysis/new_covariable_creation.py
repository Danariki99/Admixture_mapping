import sys
import subprocess
import os
import pandas as pd
import polars as pl

if __name__ == '__main__':
    # Check if the dataset argument is provided
    if len(sys.argv) != 2:
        print("Usage: python post_processing.py <dataset>")
        sys.exit(1)

    # The dataset variable is taken from the command line argument
    dataset = sys.argv[1]

    # Use the dataset variable to construct file paths
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files'
    wind_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/{dataset}/wind'
    tmp_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp'
    generic_msp_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb/ukb_hap_chr*_v2_rfmix.msp.tsv'
    old_covar_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered.phe'

    ancestry_map = {
        "AFR": 0,
        "AHG": 1,
        "EAS": 2,
        "EUR": 3,
        "NAT": 4,
        "OCE": 5,
        "SAS": 6,
        "WAS": 7
    }

    old_covar_df = pd.read_csv(old_covar_file, sep='\t')

    list_of_files = os.listdir(wind_folder)

    for wind_filename in [list_of_files[0]]:
        ancestry = wind_filename.split('_')[0]
        pheno = wind_filename.split('_')[1]
        input_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/{dataset}/ancestry_{ancestry}.vcf'
        wind_file = os.path.join(wind_folder, wind_filename)
        wind_df_1 = pd.read_csv(wind_file, sep='\t')
        chroms = wind_df_1['chr'].unique()

        for chr in [chroms[0]]:
            msp_file = generic_msp_file.replace('*', str(chr))
            output_file_1 = os.path.join(tmp_folder, 'processed_msp_1.tsv')

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

            # Read the current window file into pandas and extract chrom, start, end
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

            print(extended_wind_df)
            # Now, execute awk commands one by one for each row in extended_wind_df

            # First, output the header (second row)
            output_file_2 = os.path.join(tmp_folder, 'processed_msp_2.tsv')
            print(f"Printing header row to {output_file_2}")
            awk_header_command = f"awk 'NR == 2 {{print}}' {msp_file} > {output_file_2}"
            subprocess.run(awk_header_command, shell=True, check=True)

            # Then, run the awk command for each row of the extended window data
            for _, row in extended_wind_df.iterrows():
                # Generate the awk command for this particular row
                awk_command = f"awk '$1 == \"{row['chr']}\" && $2 == {row['start']} && $3 == {row['end']} {{print}}' {msp_file} >> {output_file_2}"
                print(f"Running awk command for row: {row['chr']} {row['start']} {row['end']}")
                subprocess.run(awk_command, shell=True, check=True)

            # Load the final processed file into a Polars DataFrame
            polars_df = pl.read_csv(output_file_2, separator=' ')

            # Convert Polars DataFrame to Pandas DataFrame
            pandas_df = polars_df.to_pandas()

            total_counts = 2 * len(pandas_df)

            # Step 1: Count occurrences of each ancestry for each row (for both haplotypes .0 and .1)
            ancestry_columns = [f'{id}.0' for id in old_covar_df['IID']] + [f'{id}.1' for id in old_covar_df['IID']]
            counts_df = pd.DataFrame(columns=ancestry_map.keys())

            for ancestry in ancestry_map.keys():
                # Count occurrences for the ancestry in both haplotypes (.0 and .1)
                counts_df[ancestry] = (
                    (pandas_df[ancestry_columns].values == ancestry_map[ancestry]).sum(axis=1)
                )

            # Step 2: Determine the ancestry with the most counts for each individual
            dominant_ancestry = counts_df.idxmax(axis=1)  # This gives the ancestry with the highest count for each row

            # Step 3: Create the ancestry indicator column and the percentage column
            pandas_df[f'{ancestry}'] = dominant_ancestry == ancestry  # Indicator: 1 if the ancestry is the dominant one
            pandas_df[f'{ancestry}_percentage'] = counts_df[ancestry] / (2 * len(pandas_df))  # Percentage of the current ancestry

            # Step 4: Merge the updated data with old_covar_df
            covar_df = pd.merge(old_covar_df, pandas_df[['IID', f'{ancestry}', f'{ancestry}_percentage']], on='IID', how='inner')

            # Step 5: Save or use the covar_df
            covar_file = os.path.join(output_folder, f'{ancestry}_{pheno}_covar_chr{chr}.tsv')
            covar_df.to_csv(covar_file, sep='\t', index=False)
