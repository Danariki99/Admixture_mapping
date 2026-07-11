import os
import sys
import subprocess
import pandas as pd
import numpy as np

# TEST version of window_covar_creation.py — same structure/logic, but runs on the
# data we pass (result_folder / data_folder) instead of the hardcoded production
# paths. Builds one per-window covariate file that adds the window LAI (local
# ancestry dosage) to the base covariates, so the conditional fine mapping can
# test the SNP (ADD) and the ancestry (LAI) in the same GLM.

if len(sys.argv) != 3:
    print("Usage: python window_covar_creation_test.py <result_folder> <data_folder>")
    sys.exit(1)

result_folder = sys.argv[1]
data_folder   = sys.argv[2]

output_folder = os.path.join(result_folder, 'wind_covar_files_new')
wind_folder   = os.path.join(result_folder, 'FUMA/wind')
# per-ancestry window-LAI VCF (each row = a window, per-sample LAI dosage)
vcf_folder    = os.path.join(result_folder, 'vcf_files')
# base covariates WITH global-ancestry proportions (from the covariate step)
base_covar_path = os.path.join(result_folder, 'covar_file/covar_proportions.phe')
tmp_folder    = os.path.join(result_folder, 'tmp')

os.makedirs(output_folder, exist_ok=True)
os.makedirs(tmp_folder, exist_ok=True)

base_covar_df = pd.read_csv(base_covar_path, sep='\t')
base_covar_df['IID'] = base_covar_df['IID'].astype(str)


def gt_to_dosage(gt):
    return gt.count('1')  # 0/0->0, 0/1->1, 1/0->1, 1/1->2


for wind_filename in os.listdir(wind_folder):
    ancestry = wind_filename.split('_')[0]
    pheno    = wind_filename.split('_')[1]
    vcf_file = os.path.join(vcf_folder, f'ancestry_{ancestry}.vcf')
    wind_df_1 = pd.read_csv(os.path.join(wind_folder, wind_filename), sep='\t')

    # sample names from VCF header (columns 10+)
    header_result = subprocess.run(
        f"grep '^#CHROM' {vcf_file}", shell=True, capture_output=True, text=True, check=True
    )
    vcf_samples = header_result.stdout.strip().split('\t')[9:]

    for chr_val in wind_df_1['chr'].unique():
        print(f"Processing {wind_filename} chr{chr_val}...")
        wind_chr = wind_df_1[wind_df_1['chr'] == chr_val]

        for _, row in wind_chr.iterrows():
            start, end = int(row['start']), int(row['end'])
            pos_str = f"{start}_{end}"

            row_result = subprocess.run(
                f"awk '!/^##/ && $1 == \"{chr_val}\" && $2 == \"{pos_str}\" {{print}}' {vcf_file}",
                shell=True, capture_output=True, text=True
            )
            if not row_result.stdout.strip():
                print(f"  Warning: window {chr_val}:{pos_str} not found in VCF, skipping")
                continue

            genotypes = row_result.stdout.strip().split('\t')[9:]
            if len(genotypes) != len(vcf_samples):
                print(f"  Warning: field count mismatch, skipping")
                continue

            lai_df = pd.DataFrame({'IID': [str(s) for s in vcf_samples],
                                   'LAI': [gt_to_dosage(g) for g in genotypes]})
            merged = pd.merge(base_covar_df, lai_df, on='IID', how='inner')

            out_path = os.path.join(
                output_folder, f'{ancestry}_{pheno}_chr{chr_val}_{start}_{end}_covar.tsv')
            merged.to_csv(out_path, sep='\t', index=False)
            print(f"  Created {out_path}")
