import os
import subprocess
import pandas as pd

vcf_file    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb.vcf.gz'
king_keep   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
wind_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
base_output = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld'
sbatch_dir  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/ukbb/SuSiE_ld'

os.makedirs(sbatch_dir, exist_ok=True)

for wind_filename in sorted(os.listdir(wind_folder)):
    if not wind_filename.endswith('_wind.txt'):
        continue

    ancestry = wind_filename.split('_')[0]
    pheno    = wind_filename.split('_')[1]

    wind_file = os.path.join(wind_folder, wind_filename)
    wind_df   = pd.read_csv(wind_file, sep='\t')

    for chr_val in wind_df['chr'].unique():
        wind_chr     = wind_df[wind_df['chr'] == chr_val]
        region_start = int(wind_chr['start'].min())
        region_end   = int(wind_chr['end'].max())

        hit_label     = f'{ancestry}_{pheno}_chr{chr_val}'
        output_folder = os.path.join(base_output, hit_label)
        os.makedirs(output_folder, exist_ok=True)

        out_prefix = os.path.join(output_folder, f'{hit_label}_ld')

        if os.path.exists(f'{out_prefix}.phased.vcor1'):
            print(f"Already done, skipping: {out_prefix}.phased.vcor1")
            continue

        print(f"{hit_label}: region {region_start}-{region_end}")

        plink_cmd = (
            f'/private/home/rsmerigl/plink2 '
            f'--vcf {vcf_file} '
            f'--chr {chr_val} '
            f'--from-bp {region_start} '
            f'--to-bp {region_end} '
            f'--keep {king_keep} '
            f'--r-phased square '
            f'--out {out_prefix}'
        )

        sbatch_file = os.path.join(sbatch_dir, f'job_{hit_label}.sh')
        out_file    = sbatch_file.replace('.sh', '.out')
        err_file    = sbatch_file.replace('.sh', '.err')

        with open(sbatch_file, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --partition=long\n')
            f.write(f'#SBATCH --job-name=susie_ld_{hit_label}\n')
            f.write(f'#SBATCH --output={out_file}\n')
            f.write(f'#SBATCH --error={err_file}\n')
            f.write('#SBATCH --nodes=1\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write('#SBATCH --time=1-00:00:00\n')
            f.write('#SBATCH --mem=128G\n')
            f.write(f'{plink_cmd}\n')

        result = subprocess.run(['sbatch', sbatch_file], capture_output=True, text=True)
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted job {job_id}: {hit_label} ({region_start}-{region_end})")
z