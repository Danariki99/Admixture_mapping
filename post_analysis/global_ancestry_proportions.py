import os
import pandas as pd
import numpy as np


ANCESTRIES = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']

counts_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/counts'
old_covar_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered.phe'
output_percentages = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/percentages.csv'
output_covar = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered_proportions.phe'


# Step 1 — Sum count files across all 22 chromosomes
samples = None
total_counts = None

for chrom in range(1, 23):
    count_file = os.path.join(counts_folder, f'final_counts_chr{chrom}.csv')
    print(f"Reading chromosome {chrom}...")

    df_chr = pd.read_csv(count_file)

    if total_counts is None:
        samples = df_chr['#IID'].astype(str)
        total_counts = df_chr[ANCESTRIES].copy().astype(np.int64)
    else:
        total_counts += df_chr[ANCESTRIES].values

print("All chromosomes loaded.")


# Step 1b — Sanity check: every sample should have the same total window count
row_sums = total_counts.sum(axis=1)
n_unique = row_sums.nunique()
print(f"Total windows per sample — min: {row_sums.min()}, max: {row_sums.max()}, unique values: {n_unique}")
if n_unique != 1:
    print(f"WARNING: not all samples have the same total window count! ({n_unique} unique values)")
else:
    print(f"OK: all samples have the same total window count ({row_sums.iloc[0]})")


# Step 2 — Compute percentages and save percentages.csv
percentages = total_counts.div(row_sums, axis=0) * 100

percentages_df = percentages.copy()
percentages_df.insert(0, 'IID', samples.values)

percentages_df.to_csv(output_percentages, index=False)
print(f"Saved percentages to {output_percentages}")


# Step 3 — Build new covariate file replacing PCs with ancestry proportions
covar_df = pd.read_csv(old_covar_file, sep='\t')

pc_cols = [c for c in covar_df.columns if c.startswith('Global_PC')]
covar_df = covar_df.drop(columns=pc_cols)

covar_df['IID'] = covar_df['IID'].astype(str)
percentages_df['IID'] = percentages_df['IID'].astype(str)

covar_df = pd.merge(covar_df, percentages_df, on='IID')

covar_df.to_csv(output_covar, sep='\t', index=False)
print(f"Saved new covariate file to {output_covar}")
