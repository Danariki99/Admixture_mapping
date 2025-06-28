import os
import subprocess
import time
import pandas as pd
import sys
import numpy as np

result_folder = sys.argv[1]

ancestry_list = ['AFR', 'EAS', 'EUR']

# Define the folder where the files are stored
hit_folder = os.path.join(result_folder, 'probabities_pipeline/models')
covar_folder = os.path.join(result_folder, 'PCA_files/PCA_covar_files')
keep_folder = os.path.join(result_folder, 'keep_files_processed')
output_folder = os.path.join(result_folder, 'probabities_pipeline/samples')
tmp_folder = os.path.join(result_folder, 'probabities_pipeline/tmp')
vcf_folder = os.path.join(result_folder, 'vcf_folder')

if not os.path.exists(tmp_folder):
    os.makedirs(tmp_folder)

hit_list = os.listdir(hit_folder)

for hit in hit_list:
    current_folder = f'{hit_folder}/{hit}'
    pheno = hit.split('_')[1]
    for ancestry in ancestry_list:
        output_folder_file = os.path.join(output_folder, f'{hit}/{ancestry}')
        os.makedirs(output_folder_file, exist_ok=True)

        keep_name = f'{ancestry}_keep'
        keep_file = f'{keep_folder}/{keep_name}.txt'
        if not os.path.exists(keep_file):
            print(f'Keep file {keep_file} does not exist, skipping.')
            continue
        keep_df = pd.read_csv(keep_file, sep='\t')
        keep_df.rename(columns={'#IID': 'IID'}, inplace=True)
        current_folder = f'{hit_folder}/{hit}/{ancestry}'
        covar_name = '_covar_'.join(['_'.join(hit.split('_')[:2]), hit.split('_')[2]])
        covar_file = f'{covar_folder}/{covar_name}_{ancestry}.tsv'
        if not os.path.exists(covar_file):
            print(f'Covariate file {covar_file} does not exist, skipping.')
            continue
        covar_df = pd.read_csv(covar_file, sep='\t')
        files_list = os.listdir(current_folder)


        extract_file = os.path.join(tmp_folder, "snps_list.txt")

        counter = 0
        with open(extract_file, "w") as f:
            for file_name in files_list:
                if file_name.endswith(".csv"):
                    snp_id = file_name.split('.')[0]
                    f.write(snp_id + "\n")
                    counter = counter+1
        if counter == 0:
            print('skipping caused by missing snps')
            continue

        print(counter)
        print("extract file: ", extract_file)
        
        # use bcftools to extract the single SNP from the chromosomve vcf file. 
        output_plink = os.path.join(tmp_folder, f'{hit}')
        chrom = hit.split('_')[-1]
        chrom_vcf_file = f'{vcf_folder}/{chrom}.vcf.gz'
        
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

            raw_df.columns = ['_'.join(col.split('_')[:-1]) if '_' in col else col for col in raw_df.columns]

            for beta_file in files_list:
                id = beta_file.split('.')[0]
                output_file = os.path.join(output_folder, f'{hit}/{ancestry}/{id}.tsv')
                beta_df = pd.read_csv(f'{current_folder}/{beta_file}')
                print(f'{current_folder}/{beta_file}')
                print(beta_df.head())
                

                interesting_columns = ['IID'] + [beta_df.columns.tolist()[-1]]

                dataset_df = covar_df[interesting_columns]
                add_df = raw_df[['IID', id]]
                dataset_df = dataset_df.merge(add_df, on='IID')
                dataset_df = dataset_df[dataset_df['IID'].isin(keep_df['IID'])]
                dataset_df.rename(columns={id: 'ADD'}, inplace=True)


                output_file = os.path.join(output_folder_file, f'{id}.tsv')

                if dataset_df.select_dtypes(include=[np.number]).applymap(np.isfinite).all().all():
                    dataset_df.to_csv(output_file, sep='\t', index=False)
                else:
                    print(f"⚠️ File non salvato (valori non numerici): {output_file}")
    
