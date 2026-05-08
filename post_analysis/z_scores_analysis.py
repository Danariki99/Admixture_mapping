import os
import pandas as pd
import numpy as np
import scipy.stats as stats


z_scores_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/z_scores/results/'
no_covar_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_no_covar_PCA/'
covar_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose/'

hits = os.listdir(z_scores_folder)

for hit in hits:
    ancestry = hit.split('_')[0]
    pheno = hit.split('_')[1]

    z_scores = pd.read_csv(os.path.join(z_scores_folder, hit,f'{hit}_output.{ancestry}.z_scores.txt'), sep='\t')

    z_scores['BETA_ABS'] = abs(np.log(z_scores['OR']))

    important_ID = z_scores.sort_values('BETA_ABS', ascending=False).iloc[0]['ID']
    print(important_ID)

    no_covar = pd.read_csv(os.path.join(no_covar_path, hit, f'{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'), sep='\t')

    covar = pd.read_csv(os.path.join(covar_path, hit, f'{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'), sep='\t')

    no_covar_row = no_covar[(no_covar['ID'] == important_ID) & (no_covar['TEST'] == 'ADD')]

    covar_row = covar[(covar['ID'] == important_ID) & (covar['TEST'] == 'ADD')]

    or_no_covar = no_covar_row['OR'].values[0]
    or_covar = covar_row['OR'].values[0]

    shift = abs(or_covar - 1) - abs(or_no_covar - 1)

    if shift < 0:
        print("the covariable is reducing significance")
    elif shift > 0:
        print("the covariable is increasing significance")
    else:
        print("the covariable is not changing significance")

    # compute the wald z_score on the betas
    beta_no_covar = np.log(no_covar_row['OR'].values[0])
    beta_covar = np.log(covar_row['OR'].values[0])

    se_no_covar = no_covar_row['LOG(OR)_SE'].values[0]
    se_covar = covar_row['LOG(OR)_SE'].values[0]

    merged_z = (beta_covar - beta_no_covar) / np.sqrt(se_no_covar**2 + se_covar**2)

    p = 2 * (1 - stats.norm.cdf(abs(merged_z)))

    print(f"the p-value of the merged z-score is {p} so the covariable is {'' if p < 0.05 else 'not'} significant")

    




