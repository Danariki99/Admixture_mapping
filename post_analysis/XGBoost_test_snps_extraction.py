import os
import subprocess
import pandas as pd

# Paths
vcf_file    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb.vcf.gz'
keep_file   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
wind_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
base_output = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/XGBoost_test_snps'
sbatch_dir  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/ukbb/XGBoost_test_snps'

os.makedirs(sbatch_dir, exist_ok=True)

for wind_filename in sorted(os.listdir(wind_folder)):
    if not wind_filename.endswith('_wind.txt'):
        continue

    ancestry = wind_filename.split('_')[0]
    pheno    = wind_filename.split('_')[1]

    wind_file = os.path.join(wind_folder, wind_filename)
    wind_df   = pd.read_csv(wind_file, sep='\t')

    for chr_val in wind_df['chr'].unique():
        wind_chr = wind_df[wind_df['chr'] == chr_val]

        region_start = int(wind_chr['start'].min())
        region_end   = int(wind_chr['end'].max())

        output_folder = os.path.join(base_output, f'{ancestry}_{pheno}_chr{chr_val}')
        os.makedirs(output_folder, exist_ok=True)

        out_prefix = os.path.join(output_folder, f'{ancestry}_{pheno}_chr{chr_val}_snps')

        if os.path.exists(f'{out_prefix}.raw'):
            print(f"Already done, skipping: {out_prefix}.raw")
            continue

        plink_cmd = (
            f'/private/home/rsmerigl/plink2 '
            f'--vcf {vcf_file} '
            f'--chr {chr_val} '
            f'--from-bp {region_start} '
            f'--to-bp {region_end} '
            f'--keep {keep_file} '
            f'--export A '
            f'--out {out_prefix}'
        )

        sbatch_file = os.path.join(sbatch_dir, f'job_{ancestry}_{pheno}_chr{chr_val}.sh')
        out_file    = sbatch_file.replace('.sh', '.out')
        err_file    = sbatch_file.replace('.sh', '.err')

        with open(sbatch_file, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --partition=long\n')
            f.write(f'#SBATCH --job-name=xgb_{ancestry}_{pheno}_{chr_val}\n')
            f.write(f'#SBATCH --output={out_file}\n')
            f.write(f'#SBATCH --error={err_file}\n')
            f.write('#SBATCH --nodes=1\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write('#SBATCH --time=1-00:00:00\n')
            f.write('#SBATCH --mem=128G\n')
            f.write(f'{plink_cmd}\n')

        result = subprocess.run(['sbatch', sbatch_file], capture_output=True, text=True)
        job_id  = result.stdout.strip().split()[-1]
        print(f"Submitted job {job_id}: {ancestry}_{pheno} chr{chr_val} ({region_start}-{region_end})")
