import os
import pandas as pd
from scipy.stats import wilcoxon

hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA/'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS']

for hit in os.listdir(hit_folder_name):
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/ratios/results/{hit}'
    os.makedirs(output_folder, exist_ok=True)
    for ancestry in ancestry_list:
        covar_in = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/ratios/covarIN_PCA/{hit}_output.{ancestry}.glm.logistic.hybrid.v2'
        covar_out = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/ratios/covarOUT_PCA/{hit}_output.{ancestry}.glm.logistic.hybrid.v2'

        df_in = pd.read_csv(covar_in, sep=' ')
        df_out = pd.read_csv(covar_out, sep=' ')

        P_in = df_in['P']
        P_out = df_out['P']

        W, p = wilcoxon(P_in, P_out)

        if p < 0.05:
            print(f'Hit: {hit}, Ancestry: {ancestry}, p-value: {p}')
