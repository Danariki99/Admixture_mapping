import os
import pandas as pd

ratios_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/ratios/results'
hits = os.listdir(ratios_path)

# Initialize an empty DataFrame to store results
out1 = pd.DataFrame(columns=[
    'Phenotype',
    'ID',
    'ancestry tested',
    'Ratio',
    'chr',
])

for hit in hits:

    ancestry = hit.split('_')[0]

    chr = hit.split('_')[-1].replace('chr', '')

    ratio_file = os.path.join(ratios_path, hit, f'{hit}_output.{ancestry}.ratios.txt')
    ratio_df = pd.read_csv(ratio_file, sep=' ')

    row_with_max_ratio = ratio_df.loc[ratio_df['P_ratio'].idxmax()]
    row_with_min_ratio = ratio_df.loc[ratio_df['P_ratio'].idxmin()]

    max_impact = max(row_with_max_ratio['P_ratio'], 1 / row_with_max_ratio['P_ratio'])
    min_impact = max(row_with_min_ratio['P_ratio'], 1 / row_with_min_ratio['P_ratio'])

    chosen_row = row_with_max_ratio if max_impact > min_impact else row_with_min_ratio

    ratio = chosen_row['P_ratio']
    snpid = chosen_row['ID']
    pheno_name = hit.split('_')[1]


    out1 = pd.concat([
        out1,
        pd.DataFrame([{
            'Phenotype': pheno_name,
            'ID': snpid,
            'ancestry tested': ancestry,  
            'Ratio': ratio,
            'chr': chr
        }])

    ], ignore_index=True)

out1.to_excel('table_ratios.xlsx', index=False)


