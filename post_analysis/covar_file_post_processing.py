import os
import pandas as pd

covar_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files'

keep_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/ukbb/keep_file.txt'

keep_df = pd.read_csv(keep_file, sep='\t')

output_folder  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_processed'

for covar_file in os.listdir(covar_folder):
    covar_df = pd.read_csv(os.path.join(covar_folder, covar_file), sep='\t')

    ancestry = covar_file.split('_')[0]

    keep_df = keep_df.rename(columns={'#IID': 'IID'})

    # Filtro il DataFrame in base agli IID presenti nel file keep
    covar_df_subset = covar_df[covar_df['IID'].isin(keep_df['IID'])]

    if covar_df_subset.isna().any().any():
        print(f'NAs present in {covar_file}')
    
    covar_df_subset.to_csv(os.path.join(output_folder, covar_file), sep='\t', index=False)