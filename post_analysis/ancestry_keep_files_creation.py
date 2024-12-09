import os
import pandas as pd
import numpy as np
import snputils as su
import time

msp_files_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb'

output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb'

# AFR=0       AHG=1   EAS=2   EUR=3   NAT=4   OCE=5   SAS=6   WAS=7

df_counts = pd.DataFrame(columns=['#IID', 'AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS'])

# Load the list of the files

msp_files = os.listdir(msp_files_folder)
msp_files = sorted(
    msp_files,
    key=lambda x: int(x.split('_chr')[1].split('_')[0])
)

for msp_file in msp_files:
    print('Processing', msp_file)

    # Load the msp file
    start_time = time.time()
    laiobj = su.MSPReader(os.path.join(msp_files_folder, msp_file)).read()

    # Initialize df_counts if empty
    if df_counts.empty:
        df_counts['#IID'] = laiobj.samples
        df_counts.iloc[:, 1:] = 0  # Initialize counts to 0

    # Update the DataFrame with counts
    n_wind, n_haplotypes = laiobj.lai.shape
    n_samp = len(laiobj.samples)

    print('start counting')

    # Vettorializzazione per evitare il ciclo annidato
    hap_to_sample = np.arange(n_haplotypes) // 2  # Map haplotype index to sample index
    ancestry_codes = laiobj.lai  # Entire LAI matrix

    ancestry_to_col = {v: i for i, v in enumerate(df_counts.columns[1:])}  # Map ancestry labels to columns

    temp_counts = np.zeros((n_samp, len(df_counts.columns) - 1), dtype=int)  # Temporary counts array

    for wind in range(n_wind):
        ancestry_labels = [laiobj.ancestry_map.get(str(code), None) for code in ancestry_codes[wind, :]]
        ancestry_indices = [ancestry_to_col.get(label, None) for label in ancestry_labels]


        valid_indices = [(hap_idx, col) for hap_idx, col in enumerate(ancestry_indices) if col is not None]
        for hap_idx, col in valid_indices:
            temp_counts[hap_to_sample[hap_idx], col] += 1

    # Add temporary counts to the DataFrame
    df_counts.iloc[:, 1:] += temp_counts

    print('Time to load the msp file:', time.time() - start_time)

# Save the final df_counts if needed
output_file = os.path.join(output_folder, 'final_counts.csv')
df_counts.to_csv(output_file, index=False)

# Generate ancestry-specific keep files
for ancestry in df_counts.columns[1:]:
    max_ancestry_samples = df_counts.loc[df_counts[ancestry] == df_counts.iloc[:, 1:].max(axis=1), '#IID']
    keep_file = os.path.join(output_folder, ancestry + '_keep.txt')
    keep_df = pd.DataFrame(max_ancestry_samples, columns=['#IID'])
    keep_df.to_csv(keep_file, index=False, header=True, sep='\t')

print('Processing complete. Final counts and keep files have been saved.')
