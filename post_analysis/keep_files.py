import os
import pandas as pd
import numpy as np

#modify this thing to take all the files in the final output folder
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/keep_files'

counts_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/counts'

counts_files = os.listdir(counts_folder)

df_counts = None

for count_file in counts_files:
    print('Processing ' + count_file)
    df = pd.read_csv(os.path.join(counts_folder, count_file))
    df = df.set_index('#IID')
    df = df.dropna()
    df = df.astype(int)
    if df_counts is None:
        df_counts = df
    else:
        df_counts = df_counts.add(df, fill_value=0)
    
df_counts = df_counts.reset_index()

# Generate ancestry-specific keep files
for ancestry in df_counts.columns[1:]:
    max_ancestry_samples = df_counts.loc[df_counts[ancestry] == df_counts.iloc[:, 1:].max(axis=1), '#IID']
    keep_file = os.path.join(output_folder, ancestry + '_keep.txt')
    keep_df = pd.DataFrame(max_ancestry_samples, columns=['#IID'])
    keep_df.to_csv(keep_file, index=False, header=True, sep='\t')



print('Processing complete. Final counts and keep files have been saved.')