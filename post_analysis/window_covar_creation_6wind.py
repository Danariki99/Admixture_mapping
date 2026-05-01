import os
import subprocess
import pandas as pd
import numpy as np


output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new_6wind'
wind_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
vcf_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/ukbb'
base_covar_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar_filtered_proportions.phe'
tmp_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp'

N_EXT = 6

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

        # Get all window positions on this chrom from ancestry VCF
        tmp_pos = os.path.join(tmp_folder, f'{ancestry}_{pheno}_{chr_val}_pos.txt')
        subprocess.run(
            f"awk '!/^##/ && $1 == \"{chr_val}\" {{print $2}}' {vcf_file} > {tmp_pos}",
            shell=True, check=True
        )

        pos_df = pd.read_csv(tmp_pos, header=None, names=['pos'])
        pos_df[['start', 'end']] = pos_df['pos'].str.split('_', expand=True)
        pos_df['start'] = pos_df['start'].astype(int)
        pos_df['end'] = pos_df['end'].astype(int)
        pos_df = pos_df.drop(columns=['pos']).sort_values('start').reset_index(drop=True)
        os.remove(tmp_pos)

        wind_chr = wind_df_1[wind_df_1['chr'] == chr_val]
        min_start = wind_chr['start'].min()
        max_end = wind_chr['end'].max()

        start_idx = pos_df[pos_df['start'] == min_start].index[0]
        end_idx = pos_df[pos_df['end'] == max_end].index[0]

        # Extend ±6 windows around the significant region
        extended = pd.concat([
            pos_df.iloc[max(0, start_idx - N_EXT):start_idx],
            pos_df[(pos_df['start'] >= min_start) & (pos_df['end'] <= max_end)],
            pos_df.iloc[end_idx + 1:end_idx + 1 + N_EXT]
        ]).drop_duplicates().reset_index(drop=True)

        # Only create covariate files for extension windows (sig ones already in wind_covar_files_new)
        sig_positions = set(zip(wind_chr['start'].astype(int), wind_chr['end'].astype(int)))
        extension_only = extended[~extended.apply(lambda r: (int(r['start']), int(r['end'])) in sig_positions, axis=1)]
        print(f"  {len(wind_chr)} sig windows (skipped) + {len(extension_only)} new extension windows")

        for _, row in extension_only.iterrows():
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
                print(f"  Warning: field count mismatch ({len(genotypes)} vs {len(vcf_samples)}), skipping")
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
