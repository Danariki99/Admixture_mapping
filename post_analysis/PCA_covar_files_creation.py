import os
import pandas as pd

ancestries = ['AFR', 'EAS', 'EUR', 'NAT', 'SAS', 'WAS']

# File paths
covar_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_processed/'
eigenvec_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_res/PCA_*.out.eigenvec'
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_covar_files'

for covar_file in os.listdir(covar_path):

    covar_df = pd.read_csv(os.path.join(covar_path, covar_file), sep='\t')

    # Load dataframes
    for ancestry in ancestries:

        eigenvec_df = pd.read_csv(eigenvec_path.replace('*', ancestry), sep='\t')

        # Rename eigenvec IID column for consistency
        eigenvec_df.rename(columns={'#IID': 'IID'}, inplace=True)

        # Keep only non-PC columns from covar_df
        covar_df_filtered = covar_df.drop(columns=[col for col in covar_df.columns if col.startswith("Global_PC")])

        covarage_sex = covar_df_filtered[['IID', 'age', 'sex', 'BMI']]
        covar_filered = covar_df_filtered.drop(columns=['age', 'sex', 'BMI'])

        # Merge covar_df with eigenvec_df on IID
        merged_df = pd.merge(pd.merge(covarage_sex, eigenvec_df, on='IID', how='inner'), covar_filered, on='IID')

        # Define the output path
        output_path = os.path.join(output_folder, f"{covar_file.split('.')[0]}_{ancestry}.tsv")

        # Save the merged dataframe
        merged_df.to_csv(output_path, sep='\t', index=False)
