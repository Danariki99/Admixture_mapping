import os
import pandas as pd
import polars as pl
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

# Just to take only one for the test
msp_files = [msp_files[-1]]

for msp_file in msp_files:
    print('Processing', msp_file)

    # Load the msp file
    start_time = time.time()
    laiobj = su.MSPReader(os.path.join(msp_files_folder, msp_file)).read()
    print('Time to load the msp file:', time.time() - start_time)

    print(len(laiobj.samples))
    print(laiobj.lai.shape)
    print(laiobj.ancestry_map)

    # Initialize df_counts if empty
    if df_counts.empty:
        df_counts['#IID'] = laiobj.samples
        df_counts.iloc[:, 1:] = 0  # Initialize counts to 0

    # Update the DataFrame with counts
    n_wind, n_haplotypes = laiobj.lai.shape
    n_samp = len(laiobj.samples)

    for wind in range(n_wind):
        for hap_idx in range(n_haplotypes):
            sample_idx = hap_idx // 2  # Derive sample index from haplotype index
            sample = laiobj.samples[sample_idx]

            ancestry_code = laiobj.lai[wind, hap_idx]
            ancestry_label = laiobj.ancestry_map.get(ancestry_code, None)

            if ancestry_label in df_counts.columns:
                df_counts.loc[df_counts['#IID'] == sample, ancestry_label] += 1

# Save the final df_counts if needed
output_file = os.path.join(output_folder, 'final_counts.csv')
df_counts.to_csv(output_file, index=False)
'''
# Generate ancestry-specific keep files
for ancestry in df_counts.columns[1:]:
    max_ancestry_samples = df_counts.loc[df_counts[ancestry] == df_counts.iloc[:, 1:].max(axis=1), '#IID']
    keep_file = os.path.join(output_folder, ancestry + '_keep.txt')
    keep_df = pd.DataFrame(max_ancestry_samples, columns=['#IID'])
    keep_df.to_csv(keep_file, index=False, header=True, sep='\t')

print('Processing complete. Final counts and keep files have been saved.')
'''