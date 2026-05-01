import os
import pandas as pd

keep_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'

keep_folder   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files_old/ukbb/keep_files'
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/keep_files_processed'

os.makedirs(output_folder, exist_ok=True)

keep_df = pd.read_csv(keep_file, sep='\t')

for anc_file in os.listdir(keep_folder):
    anc_df = pd.read_csv(os.path.join(keep_folder, anc_file), sep='\t')
    filtered = anc_df[anc_df['#IID'].isin(keep_df['#IID'])]
    filtered.to_csv(os.path.join(output_folder, anc_file), sep='\t', index=False)
    print(f'{anc_file}: {len(anc_df)} → {len(filtered)} samples')
