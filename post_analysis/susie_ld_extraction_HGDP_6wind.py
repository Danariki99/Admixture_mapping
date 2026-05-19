import os
import subprocess
import pandas as pd
from collections import defaultdict

# 1380-sample reference panel (hg19, ADMIXTURE ≥95% pure individuals)
hgdp_bfile             = '/private/groups/ioannidislab/galangal_dirs/ref_1kg_hgdp_sgdp/beagle_1kg_hgdp_sgdp_ref_panel_hg19_pure'
rfmix_map              = '/private/home/rsmerigl/codes/cleaned_codes/usefull_panel_data/rfmix_sample_map.tsv'
fine_mapping_folder    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
fine_mapping_6w_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new_6wind'
covar_folder_6wind     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new_6wind'
base_output            = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld_HGDP_6wind'
sbatch_dir             = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/ukbb/SuSiE_ld_HGDP_6wind'

os.makedirs(sbatch_dir, exist_ok=True)
os.makedirs(base_output, exist_ok=True)

# Build per-ancestry keep files from the 1380-sample panel map
meta = pd.read_csv(rfmix_map, sep='\t').rename(columns={'#Sample': 'sample_id'})
keep_files = {}
for anc, grp in meta.groupby('Panel'):
    keep_path = os.path.join(base_output, f'keep_{anc}_HGDP.txt')
    with open(keep_path, 'w') as f:
        f.write('#IID\n')
        for sid in sorted(grp['sample_id'].astype(str)):
            f.write(f'{sid}\n')
    keep_files[anc] = keep_path
    print(f"Keep file: {anc} — {len(grp)} samples")


def load_zscores_positions(hit_label):
    """Load z-score positions (hg19) from both fine-mapping folders."""
    rows = []
    for folder in [fine_mapping_folder, fine_mapping_6w_folder]:
        hit_path = os.path.join(folder, hit_label)
        if not os.path.isdir(hit_path):
            continue
        for filename in sorted(os.listdir(hit_path)):
            if not filename.endswith('.glm.logistic.hybrid'):
                continue
            filepath = os.path.join(hit_path, filename)
            try:
                df = pd.read_csv(filepath, sep='\t', usecols=[0, 1, 10])
                df.columns = ['CHROM', 'POS', 'TEST']
                rows.append(df[df['TEST'] == 'ADD'][['CHROM', 'POS']])
            except Exception:
                continue
    if not rows:
        return pd.DataFrame(columns=['CHROM', 'POS'])
    return pd.concat(rows).drop_duplicates().reset_index(drop=True)


# Group covar files by hit label to get the extended region boundaries
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
    chr_part  = parts[2]   # e.g. 'chr5'
    start     = int(parts[3])
    end       = int(parts[4])
    chr_val   = chr_part.replace('chr', '')
    hit_label = f'{ancestry}_{pheno}_{chr_part}'

    hit_regions[hit_label]['starts'].append(start)
    hit_regions[hit_label]['ends'].append(end)
    hit_regions[hit_label]['ancestry'] = ancestry
    hit_regions[hit_label]['chr_val']  = chr_val

# Submit one job per hit
for hit_label, info in sorted(hit_regions.items()):
    ancestry     = info['ancestry']
    chr_val      = info['chr_val']
    region_start = min(info['starts'])
    region_end   = max(info['ends'])

    if ancestry not in keep_files:
        print(f"No keep file for ancestry {ancestry}, skipping")
        continue

    output_folder = os.path.join(base_output, hit_label)
    os.makedirs(output_folder, exist_ok=True)

    out_prefix = os.path.join(output_folder, f'{hit_label}_ld_HGDP_6wind')

    if os.path.exists(f'{out_prefix}.phased.vcor1'):
        print(f"Already done, skipping: {out_prefix}.phased.vcor1")
        continue

    pos_df = load_zscores_positions(hit_label)
    if pos_df.empty:
        print(f"No z-score positions for {hit_label}, skipping")
        continue

    # Write extract file with hg19 CHR:POS IDs — no liftover needed
    extract_path = os.path.join(output_folder, 'extract_positions.txt')
    with open(extract_path, 'w') as f:
        for _, row in pos_df.iterrows():
            f.write(f"{row['CHROM']}:{int(row['POS'])}\n")

    print(f"{hit_label}: {len(pos_df)} positions, region hg19 {region_start}-{region_end}")

    plink_cmd = (
        f'/private/home/rsmerigl/plink2 '
        f'--bfile {hgdp_bfile} '
        f'--chr {chr_val} '
        f'--from-bp {region_start} '
        f'--to-bp {region_end} '
        f'--keep {keep_files[ancestry]} '
        f'--set-all-var-ids @:# '
        f'--extract {extract_path} '
        f'--force-intersect '
        f'--r-phased square '
        f'--out {out_prefix}'
    )

    sbatch_file = os.path.join(sbatch_dir, f'job_{hit_label}_HGDP_6wind.sh')
    out_file    = sbatch_file.replace('.sh', '.out')
    err_file    = sbatch_file.replace('.sh', '.err')

    with open(sbatch_file, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --partition=long\n')
        f.write(f'#SBATCH --job-name=susie_hgdp6_{hit_label}\n')
        f.write(f'#SBATCH --output={out_file}\n')
        f.write(f'#SBATCH --error={err_file}\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write('#SBATCH --time=4:00:00\n')
        f.write('#SBATCH --mem=16G\n')
        f.write(f'{plink_cmd}\n')

    result = subprocess.run(['sbatch', sbatch_file], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1]
    print(f"Submitted job {job_id}: {hit_label} (hg19 {region_start}-{region_end})")
