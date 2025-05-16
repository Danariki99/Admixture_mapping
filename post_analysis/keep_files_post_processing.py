import os
import pandas as pd


keep_file_old = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/ukbb/keep_file.txt'

keep_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/keep_files'

output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/keep_files_processed'

old_keep_df = pd.read_csv(keep_file_old, sep='\t')


for keep_file in os.listdir(keep_folder):
    keep_df = pd.read_csv(os.path.join(keep_folder, keep_file), sep='\t')

    keep_file_subset = keep_df[keep_df['#IID'].isin(old_keep_df['#IID'])]

    keep_file_subset.to_csv(os.path.join(output_folder, keep_file), sep='\t', index=False)
