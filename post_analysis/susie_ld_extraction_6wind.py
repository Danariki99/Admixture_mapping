import os
import subprocess
import pandas as pd
from collections import defaultdict

vcf_file           = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb.vcf.gz'
king_keep          = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
covar_folder_6wind = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new_6wind'
base_output        = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld_6wind'
sbatch_dir         = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/ukbb/SuSiE_ld_6wind'

os.makedirs(sbatch_dir, exist_ok=True)
os.makedirs(base_output, exist_ok=True)

hit_regions = defaultdict(lambda: {'starts': [], 'ends': [], 'ancestry': None, 'chr_val': None})

for fname in sorted(os.listdir(covar_folder_6wind)):
    if not fname.endswith('_covar.tsv'):
        continue
    base  = fname.replace('_covar.tsv', '')
    parts = base.split('_')
    if len(parts) < 5:
        continue
    ancestry  = parts[0]
    pheno     = parts[1]
    chr_part  = parts[2]
    start     = int(parts[3])
    end       = int(parts[4])
    chr_val   = chr_part.replace('chr', '')
    hit_label = f'{ancestry}_{pheno}_{chr_part}'

    hit_regions[hit_label]['starts'].append(start)
    hit_regions[hit_label]['ends'].append(end)
    hit_regions[hit_label]['ancestry'] = ancestry
    hit_regions[hit_label]['chr_val']  = chr_val

for hit_label, info in sorted(hit_regions.items()):
    chr_val      = info['chr_val']
    region_start = min(info['starts'])
    region_end   = max(info['ends'])

    output_folder = os.path.join(base_output, hit_label)
    os.makedirs(output_folder, exist_ok=True)

    out_prefix = os.path.join(output_folder, f'{hit_label}_ld_6wind')

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

    sbatch_file = os.path.join(sbatch_dir, f'job_{hit_label}_6wind.sh')
    out_file    = sbatch_file.replace('.sh', '.out')
    err_file    = sbatch_file.replace('.sh', '.err')

    with open(sbatch_file, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --partition=long\n')
        f.write(f'#SBATCH --job-name=susie_ld6_{hit_label}\n')
        f.write(f'#SBATCH --output={out_file}\n')
        f.write(f'#SBATCH --error={err_file}\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write('#SBATCH --time=1-00:00:00\n')
        f.write('#SBATCH --mem=128G\n')
        f.write(f'{plink_cmd}\n')

    result = subprocess.run(['sbatch', sbatch_file], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1]
    print(f"Submitted job {job_id}: {hit_label}")
