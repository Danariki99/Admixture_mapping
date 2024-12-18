import os
import pandas as pd

covar_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files'

keep_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/keep_files_processed'

output_folder  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_processed'

for covar_file in os.listdir(covar_folder):
    covar_df = pd.read_csv(os.path.join(covar_folder, covar_file), sep='\t')

    ancestry = covar_file.split('_')[0]

    keep_file = os.path.join(keep_folder, f'{ancestry}_keep.txt')

    keep_df = pd.read_csv(keep_file, sep='\t')

    keep_df = keep_df.rename(columns={'#IID': 'IID'})

    # Filtro il DataFrame in base agli IID presenti nel file keep
    covar_df_subset = covar_df[covar_df['IID'].isin(keep_df['IID'])]

    covar_df_subset.to_csv(os.path.join(output_folder, covar_file), sep='\t', index=False)