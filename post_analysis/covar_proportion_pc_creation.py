import pandas as pd

proportions_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered_proportions.phe'
old_covar_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered.phe'
output_covar = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered_proportions_pc.phe'

proportions_df = pd.read_csv(proportions_file, sep='\t')
proportions_df['IID'] = proportions_df['IID'].astype(str)

old_covar_df = pd.read_csv(old_covar_file, sep='\t')
old_covar_df['IID'] = old_covar_df['IID'].astype(str)
pc3_to_10 = old_covar_df[['IID'] + [f'Global_PC{i}' for i in range(3, 11)]]

result = pd.merge(proportions_df, pc3_to_10, on='IID')

result.to_csv(output_covar, sep='\t', index=False)
print(f"Saved to {output_covar}")
print(f"Shape: {result.shape}")
print(f"Columns: {list(result.columns)}")
