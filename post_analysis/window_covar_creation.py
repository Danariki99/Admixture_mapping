import os
import subprocess
import pandas as pd
import numpy as np


output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new'
wind_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
vcf_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/ukbb'
base_covar_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered_proportions.phe'
tmp_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp'

os.makedirs(output_folder, exist_ok=True)
os.makedirs(tmp_folder, exist_ok=True)

base_covar_df = pd.read_csv(base_covar_path, sep='\t')
base_covar_df['IID'] = base_covar_df['IID'].astype(str)

def gt_to_dosage(gt):
    return gt.count('1')  # 0/0->0, 0/1->1, 1/0->1, 1/1->2

for wind_filename in os.listdir(wind_folder):
    ancestry = wind_filename.split('_')[0]
    pheno = wind_filename.split('_')[1]
    vcf_file = os.path.join(vcf_folder, f'ancestry_{ancestry}.vcf')
    wind_file = os.path.join(wind_folder, wind_filename)
    wind_df_1 = pd.read_csv(wind_file, sep='\t')

    # Get sample names from VCF header (columns 10+)
    header_result = subprocess.run(
        f"grep '^#CHROM' {vcf_file}",
        shell=True, capture_output=True, text=True, check=True
    )
    vcf_samples = header_result.stdout.strip().split('\t')[9:]

    for chr_val in wind_df_1['chr'].unique():
        print(f"Processing {wind_filename} chr{chr_val}...")

        # Only use the significant windows from the FUMA wind file
        wind_chr = wind_df_1[wind_df_1['chr'] == chr_val]

        for _, row in wind_chr.iterrows():
            start = int(row['start'])
            end = int(row['end'])
            pos_str = f"{start}_{end}"

            # Extract the full LAI row for this window
            row_result = subprocess.run(
                f"awk '!/^##/ && $1 == \"{chr_val}\" && $2 == \"{pos_str}\" {{print}}' {vcf_file}",
                shell=True, capture_output=True, text=True
            )

            if not row_result.stdout.strip():
                print(f"  Warning: window {chr_val}:{pos_str} not found in VCF, skipping")
                continue

            fields = row_result.stdout.strip().split('\t')
            genotypes = fields[9:]

            if len(genotypes) != len(vcf_samples):
                print(f"  Warning: field count mismatch ({len(genotypes)} genotypes vs {len(vcf_samples)} samples), skipping")
                continue

            dosages = [gt_to_dosage(g) for g in genotypes]

            lai_df = pd.DataFrame({'IID': vcf_samples, 'LAI': dosages})
            lai_df['IID'] = lai_df['IID'].astype(str)

            merged = pd.merge(base_covar_df, lai_df, on='IID', how='inner')

            out_path = os.path.join(
                output_folder,
                f'{ancestry}_{pheno}_chr{chr_val}_{start}_{end}_covar.tsv'
            )
            merged.to_csv(out_path, sep='\t', index=False)
            print(f"  Created {out_path}")
