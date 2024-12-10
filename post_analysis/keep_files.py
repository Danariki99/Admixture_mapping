import os
import pandas as pd
import numpy as np

#modify this thing to take all the files in the final output folder
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb'

df_counts = pd.DataFrame(columns=['#IID', 'AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS'])

# Generate ancestry-specific keep files
for ancestry in df_counts.columns[1:]:
    max_ancestry_samples = df_counts.loc[df_counts[ancestry] == df_counts.iloc[:, 1:].max(axis=1), '#IID']
    keep_file = os.path.join(output_folder, ancestry + '_keep.txt')
    keep_df = pd.DataFrame(max_ancestry_samples, columns=['#IID'])
    keep_df.to_csv(keep_file, index=False, header=True, sep='\t')

print('Processing complete. Final counts and keep files have been saved.')