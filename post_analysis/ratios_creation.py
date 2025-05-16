import os
import subprocess

R_file = "../sonia_scripts/script.merge.R"

hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA/'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS']

for hit in os.listdir(hit_folder_name):
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/ratios/results/{hit}'
    os.makedirs(output_folder, exist_ok=True)
    for ancestry in ancestry_list:
        covar_in = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/ratios/covarIN_PCA/{hit}_output.{ancestry}.glm.logistic.hybrid.v2'
        covar_out = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/ratios/covarOUT_PCA/{hit}_output.{ancestry}.glm.logistic.hybrid.v2'
        output_file = f'{output_folder}/{hit}_output.{ancestry}.ratios.txt'
        subprocess.run(['Rscript', R_file, covar_out, covar_in, output_file])
        print(f'Finished {hit} {ancestry}')