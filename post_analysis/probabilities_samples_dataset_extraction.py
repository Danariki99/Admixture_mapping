import os
import subprocess
import time
import pandas as pd

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS', 'NAT']

# Define the folder where the files are stored
hit_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/models'
covar_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_covar_files'
keep_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/keep_files_processed/'
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/samples'
tmp_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/tmp'
vcf_folder = '/private/groups/ioannidislab/ukbb24983/hap/pgen'

hit_list = os.listdir(hit_folder)

for hit in hit_list:
    current_folder = f'{hit_folder}/{hit}'
    pheno = hit.split('_')[1]
    for ancestry in ancestry_list:
        output_folder_file = os.path.join(output_folder, f'{hit}/{ancestry}')
        os.makedirs(output_folder_file, exist_ok=True)

        keep_name = f'{ancestry}_keep'
        keep_file = f'{keep_folder}/{keep_name}.txt'
        keep_df = pd.read_csv(keep_file, sep='\t')
        keep_df.rename(columns={'#IID': 'IID'}, inplace=True)
        current_folder = f'{hit_folder}/{hit}/{ancestry}'
        covar_name = '_covar_'.join(['_'.join(hit.split('_')[:2]), hit.split('_')[2]])
        covar_file = f'{covar_folder}/{covar_name}_{ancestry}.tsv'
        covar_df = pd.read_csv(covar_file, sep='\t')
        files_list = os.listdir(current_folder)


        extract_file = os.path.join(tmp_folder, "snps_list.txt")

        with open(extract_file, "w") as f:
            for file_name in files_list:
                if file_name.endswith(".csv"):
                    snp_id = file_name.split('.')[0]
                    f.write(snp_id + "\n")
        
        # use bcftools to extract the single SNP from the chromosomve vcf file. 
        output_plink = os.path.join(tmp_folder, f'{hit}')
        chrom = hit.split('_')[-1]
        chrom_vcf_file = f'{vcf_folder}/ukb_hap_{chrom}_v2.vcf.gz'
        
        # Comando Plink
        plink_command = [
            "/private/home/rsmerigl/plink2",
            "--vcf", chrom_vcf_file,
            "--extract", extract_file,
            "--keep", keep_file,
            "--make-bed",
            "--out", output_plink
        ]
        plink_numeric_command = [
            "/private/home/rsmerigl/plink2",
            "--bfile", output_plink,  # Prefisso generato dal comando precedente
            "--recode", "A",
            "--out", os.path.join(tmp_folder, "genotypes_numeric")
        ]
        if len(files_list) != 0:
            subprocess.run(plink_command, check=True)
            subprocess.run(plink_numeric_command, check=True)

            raw_df = pd.read_csv(os.path.join(tmp_folder, "genotypes_numeric.raw"), sep='\t')

            raw_df.columns = [col.split('_')[0] if '_' in col else col for col in raw_df.columns]

            for beta_file in files_list:
                id = beta_file.split('.')[0]
                output_file = os.path.join(output_folder, f'{hit}/{ancestry}/{id}.tsv')
                beta_df = pd.read_csv(f'{current_folder}/{beta_file}')

                interesting_columns = ['IID'] + beta_df.columns.tolist()[4:]

                dataset_df = covar_df[interesting_columns]

                add_df = raw_df[['IID', id]]
                dataset_df = dataset_df.merge(add_df, on='IID')
                dataset_df = dataset_df[dataset_df['IID'].isin(keep_df['IID'])]
                dataset_df.rename(columns={id: 'ADD'}, inplace=True)


                output_file = os.path.join(output_folder_file, f'{id}.tsv')
                dataset_df.to_csv(output_file, sep='\t', index=False)
    
