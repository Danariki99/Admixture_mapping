# to add as part of the supplementary information, we should include a table (excel)
# with mean age + sd in cases and controls per phenotypes, as well %women.
# I can help if you need!

import os
import pandas as pd
import openpyxl
import requests
import numpy as np

# Paths to important files and directories
significant_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_no_covar_PCA/'

# Initialize an empty DataFrame to store results
out1 = pd.DataFrame(columns=[
    'Phenotype',
    'ID',
    'ancestry tested',
    'OR (CI = 95%)',
    'p value',
    'chr',
])

df_first_batch = pd.read_excel("ukbb_v1.xlsx", sheet_name="first_batch")

folder_list = os.listdir(significant_path)

for fold in folder_list:
    ancestry = fold.split('/')[-1].split('_')[0]
    pheno = fold.split('/')[-1].split('_')[1]
    chrom = fold.split('/')[-1].split('_')[2]

    pheno_name = df_first_batch[df_first_batch['ID'] == pheno]['ID2']

    result_path = os.path.join(significant_path, fold, f'{fold}_output.{fold.split("_")[0]}.{pheno}.glm.logistic.hybrid')
    chr = int(chrom.replace('chr', ''))

    df = pd.read_table(result_path, sep='\t')

    filtered_df = df[df['#CHROM'] == chr]
    index_min = filtered_df['P'].idxmin()
    row_with_min_p = filtered_df.loc[index_min]
    p = row_with_min_p['P']
    oddr = row_with_min_p['OR']
    snpid = row_with_min_p['ID']

    out1 = pd.concat([
        out1,
        pd.DataFrame([{
            'Phenotype': pheno_name.iloc[0] if not pheno_name.empty else 'Unknown',  # Nome del fenotipo
            'ID': snpid,
            'ancestry tested': ancestry,  
            'OR (CI = 95%)': f'{oddr} ({row_with_min_p["L95"]}, {row_with_min_p["U95"]})', 
            'p value': p,
            'chr': chr
        }])

    ], ignore_index=True)

out1.to_excel('table_hits_snps_no_covar.xlsx', index=False)
