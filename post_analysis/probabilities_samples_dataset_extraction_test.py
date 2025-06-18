import os
import subprocess
import pandas as pd
import sys

result_folder = sys.argv[1]

# Define ancestry list
ancestry_list = ['AFR', 'EAS', 'EUR']

# Define the folder where the model files are stored
hit_folder = os.path.join(result_folder, 'fine_mapping_verbose')

# Define paths to other necessary folders
covar_folder = os.path.join(result_folder, 'PCA_files', 'PCA_covar_files')
keep_folder = os.path.join(result_folder, 'keep_files_processed')
output_folder = os.path.join(result_folder, 'fine_mapping_samples')
tmp_folder = os.path.join(result_folder, 'tmp')
vcf_folder = os.path.join(result_folder, 'vcf_folder')

# Ensure output and temporary directories exist
os.makedirs(output_folder, exist_ok=True)
os.makedirs(tmp_folder, exist_ok=True)

# Iterate over each hit in the hit folder
for hit in os.listdir(hit_folder):
    current_folder = os.path.join(hit_folder, hit)
    pheno = hit.split('_')[1]
    chrom = hit.split('_')[-1].replace('chr', '')

    for ancestry in ancestry_list:
        output_folder_file = os.path.join(output_folder, hit, ancestry)
        os.makedirs(output_folder_file, exist_ok=True)

        keep_file = os.path.join(keep_folder, f'{ancestry}_keep.txt')
        keep_df = pd.read_csv(keep_file, sep='\t')
        keep_df.rename(columns={'#IID': 'IID'}, inplace=True)

        covar_file = os.path.join(covar_folder, f'covar_{pheno}_{ancestry}.tsv')
        covar_df = pd.read_csv(covar_file, sep='\t')

        current_ancestry_folder = os.path.join(current_folder, ancestry)
        files_list = [f for f in os.listdir(current_ancestry_folder) if f.endswith('.csv')]

        if not files_list:
            continue

        extract_file = os.path.join(tmp_folder, "snps_list.txt")
        with open(extract_file, "w") as f:
            for file_name in files_list:
                snp_id = file_name.split('.')[0]
                f.write(snp_id + "\n")

        output_plink = os.path.join(tmp_folder, f'{hit}')
        chrom_vcf_file = os.path.join(vcf_folder, f'chr{chrom}.vcf.gz')

        # Run PLINK to extract SNPs
        plink_command = [
            "/private/home/rsmerigl/plink2",
            "--vcf", chrom_vcf_file,
            "--extract", extract_file,
            "--keep", keep_file,
            "--make-bed",
            "--out", output_plink
        ]
        subprocess.run(plink_command, check=True)

        # Convert to numeric format
        plink_numeric_command = [
            "/private/home/rsmerigl/plink2",
            "--bfile", output_plink,
            "--recode", "A",
            "--out", os.path.join(tmp_folder, "genotypes_numeric")
        ]
        subprocess.run(plink_numeric_command, check=True)

        raw_df = pd.read_csv(os.path.join(tmp_folder, "genotypes_numeric.raw"), sep='\t')
        raw_df.columns = [col.split('_')[0] if '_' in col else col for col in raw_df.columns]

        for beta_file in files_list:
            snp_id = beta_file.split('.')[0]
            beta_df = pd.read_csv(os.path.join(current_ancestry_folder, beta_file))

            interesting_columns = ['IID'] + beta_df.columns.tolist()[4:]
            dataset_df = covar_df[interesting_columns]

            add_df = raw_df[['IID', snp_id]]
            dataset_df = dataset_df.merge(add_df, on='IID')
            dataset_df = dataset_df[dataset_df['IID'].isin(keep_df['IID'])]
            dataset_df.rename(columns={snp_id: 'ADD'}, inplace=True)

            output_file = os.path.join(output_folder_file, f'{snp_id}.tsv')
            dataset_df.to_csv(output_file, sep='\t', index=False)
