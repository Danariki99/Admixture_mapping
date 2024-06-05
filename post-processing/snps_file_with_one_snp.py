import os
import pandas as pd
import numpy as np

# Create a function to find the closest SNP to the middle of a given window
def find_closest_snp(window):
    snps_in_window = snps_df[snps_df['pos'].between(window['POS'], window['end_POS'])]
    if not snps_in_window.empty:
        closest_snp = snps_in_window.loc[(snps_in_window['pos']-window['mid']).abs().idxmin()]
        return pd.Series([closest_snp['CHR'], closest_snp['pos'], closest_snp['rfid'], window['P'],  window['Phenotype'], window['Ancestry'], window['POS'], window['end_POS']])
    return pd.Series([np.nan, None, None, None, None, None])

type_list = ['no_BMI', 'BMI']

data_folder = '/private/groups/ioannidislab/smeriglio/output/*/'
for type in type_list:
    # Load the snps file
    snps_file = os.path.join(data_folder.replace('*', type), "significant_SNPs.txt")
    snps_df = pd.read_csv(snps_file, sep='\t')
    snps_df.rename(columns={'chr_name': 'CHR'}, inplace=True)
    snps_df.rename(columns={'refsnp_id': 'rfid'}, inplace=True)
    snps_df.rename(columns={'chrom_start': 'pos'}, inplace=True)
    snps_df = snps_df.drop(columns=['allele'])

    window_file = os.path.join(data_folder.replace('*', type), "significant_positions.tsv")
    window_df = pd.read_csv(window_file, sep='\t')

    # Convert the windows to intervals and calculate the mid point
    window_df['window'] = pd.IntervalIndex.from_arrays(window_df['POS'], window_df['end_POS'], closed='both')
    window_df['mid'] = window_df['window'].map(lambda x: x.mid)

    # Apply the function to each window and create a new dataframe
    snp_info_df = window_df.apply(find_closest_snp, axis=1)
    snp_info_df.columns = ['chr', 'pos', 'rfid', 'P', 'phenotype', 'ancestry', 'start', 'end']

    # Save the new file
    output_file = os.path.join(data_folder.replace('*', type), "significant_SNPs_with_single_P_values.txt")
    snp_info_df.to_csv(output_file, sep='\t', index=False)