import pandas as pd
import os

hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS']

hits_list = os.listdir(hit_folder_name)

for hit in hits_list:
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/z_scores/results/{hit}'
    os.makedirs(output_folder, exist_ok=True)
    pheno = hit.split('_')[1]
    for ancestry in ancestry_list:
        file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'

        df = pd.read_csv(file, sep='\t')

        df = df[df['TEST'] == hit.split('_')[0]]

        df['Z_STAT'] = df['Z_STAT'].abs()

        df.sort_values(by='Z_STAT', ascending=False, inplace=True)

        df.to_csv(f'{output_folder}/{hit}_output.{ancestry}.z_scores.txt', sep='\t', index=False)


