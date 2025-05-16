import pandas as pd
import os

hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS']

hits_list = os.listdir(hit_folder_name)

for hit in hits_list:
    pheno = hit.split('_')[1]
    for ancestry in ancestry_list:
        file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'

        df = pd.read_csv(file, sep='\t')

        df_add = df[df['TEST'] == 'ADD']

        df_ancestry = df[df['TEST'] == hit.split('_')[0]]

        significant_snps = df_add[df_add['P'] < 0.05]

        if len(significant_snps) == 0:
            print('No significant SNPs for', hit, ancestry, '\n')
        else:
            ids = significant_snps['ID'].tolist()

            mean_z = df_ancestry[df_ancestry['ID'].isin(ids)]['Z_STAT'].mean()
            std_z = df_ancestry[df_ancestry['ID'].isin(ids)]['Z_STAT'].std()

            print(hit, ancestry, mean_z, std_z, '\n')


