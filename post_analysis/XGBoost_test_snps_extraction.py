import os
import subprocess
import pandas as pd

# Paths
vcf_file    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb.vcf.gz'
keep_file   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
wind_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
base_output = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/XGBoost_test_snps'

# EAS hit
wind_filename = 'EAS_HC1158_wind.txt'
ancestry = wind_filename.split('_')[0]   # EAS
pheno    = wind_filename.split('_')[1]   # HC1158

wind_file = os.path.join(wind_folder, wind_filename)
wind_df   = pd.read_csv(wind_file, sep='\t')

for chr_val in wind_df['chr'].unique():
    wind_chr = wind_df[wind_df['chr'] == chr_val]

    output_folder = os.path.join(base_output, f'{ancestry}_{pheno}_chr{chr_val}')
    os.makedirs(output_folder, exist_ok=True)

    for _, row in wind_chr.iterrows():
        window_start = int(row['start'])
        window_end   = int(row['end'])

        out_prefix = os.path.join(output_folder, f'{ancestry}_{pheno}_chr{chr_val}_{window_start}_{window_end}_snps')

        print(f"Extracting SNPs for {ancestry}_{pheno} chr{chr_val}:{window_start}-{window_end}...")

        plink_command = [
            '/private/home/rsmerigl/plink2',
            '--vcf',      vcf_file,
            '--chr',      str(chr_val),
            '--from-bp',  str(window_start),
            '--to-bp',    str(window_end),
            '--keep',     keep_file,
            '--export',   'A',
            '--out',      out_prefix,
        ]

        subprocess.run(plink_command, check=True)
        print(f"  Done. Output: {out_prefix}.raw")
