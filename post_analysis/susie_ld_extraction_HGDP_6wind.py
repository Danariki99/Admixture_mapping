import os
import subprocess
import tempfile
import pandas as pd
from collections import defaultdict

hgdp_vcf_pattern      = '/private/home/cwshanks/directory/data/hgdp_1kgp/hgdp1kgp_chr{chr}.filtered.SNV_INDEL.phased.shapeit5.vcf.gz'
hgdp_meta             = '/private/home/cwshanks/directory/data/hgdp_1kgp/final_sample_populations_extended.tsv'
fine_mapping_folder   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
fine_mapping_6w_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new_6wind'
covar_folder_6wind    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new_6wind'
base_output           = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld_HGDP_6wind'
sbatch_dir            = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/sbatch_files/ukbb/SuSiE_ld_HGDP_6wind'
chain_file            = '/private/groups/ioannidislab/cwshanks/reference_genomes/hg19ToHg38.over.chain.gz'

ANCESTRY_MAP = {
    'AFR': 'AFR', 'EAS': 'EAS', 'SAS': 'CSA', 'WAS': 'MID',
    'AHG': 'CSA', 'NAT': 'AMR', 'OCE': 'OCE', 'EUR': 'EUR',
}

os.makedirs(sbatch_dir, exist_ok=True)
os.makedirs(base_output, exist_ok=True)

# Build per-ancestry keep files from HGDP metadata
meta = pd.read_csv(hgdp_meta, sep='\t')
keep_files = {}
for your_anc, hgdp_region in ANCESTRY_MAP.items():
    samples = meta[meta['genetic_region'] == hgdp_region]['sample_id'].astype(str)
    if len(samples) == 0:
        continue
    keep_path = os.path.join(base_output, f'keep_{your_anc}_HGDP.txt')
    with open(keep_path, 'w') as f:
        f.write('#IID\n')
        for sid in sorted(samples):
            f.write(f'{sid}\n')
    keep_files[your_anc] = keep_path
    print(f"Keep file: {your_anc} ({hgdp_region}) — {len(samples)} samples")


def load_zscores_positions(hit_label):
    """Load z-score positions from both fine-mapping folders (orig + 6wind extension)."""
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


def liftover_positions(pos_df, output_folder):
    bed_in  = os.path.join(output_folder, 'positions_hg19.bed')
    bed_out = os.path.join(output_folder, 'positions_hg38.bed')

    with open(bed_in, 'w') as f:
        for _, row in pos_df.iterrows():
            chrom = str(row['CHROM'])
            pos   = int(row['POS'])
            f.write(f"{chrom}\t{pos-1}\t{pos}\t{chrom}:{pos}\n")

    result = subprocess.run(
        ['CrossMap', 'bed', chain_file, bed_in, bed_out],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  CrossMap error: {result.stderr[:200]}")
        return []

    hg38_ids = []
    if os.path.exists(bed_out):
        with open(bed_out) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    hg38_ids.append(f"{parts[0]}:{int(parts[2])}")

    extract_path = os.path.join(output_folder, 'extract_positions.txt')
    with open(extract_path, 'w') as f:
        for snp_id in hg38_ids:
            f.write(f"{snp_id}\n")

    return hg38_ids


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
        print(f"No HGDP keep file for ancestry {ancestry}, skipping")
        continue

    output_folder = os.path.join(base_output, hit_label)
    os.makedirs(output_folder, exist_ok=True)

    out_prefix = os.path.join(output_folder, f'{hit_label}_ld_HGDP_6wind')

    if os.path.exists(f'{out_prefix}.phased.vcor1'):
        print(f"Already done, skipping: {out_prefix}.phased.vcor1")
        continue

    vcf_file = hgdp_vcf_pattern.format(chr=chr_val)
    if not os.path.exists(vcf_file):
        print(f"VCF not found for chr{chr_val}: {vcf_file}")
        continue

    pos_df = load_zscores_positions(hit_label)
    if pos_df.empty:
        print(f"No z-score positions for {hit_label}, skipping")
        continue

    print(f"{hit_label}: lifting {len(pos_df)} positions hg19→hg38 (region {region_start}-{region_end})...")
    hg38_ids = liftover_positions(pos_df, output_folder)

    if not hg38_ids:
        print(f"  No positions lifted, skipping")
        continue

    print(f"  {len(hg38_ids)}/{len(pos_df)} positions lifted successfully")

    # Liftover region boundaries
    region_bed = os.path.join(output_folder, 'region_hg38.bed')
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
        tmp.write(f"{chr_val}\t{region_start-1}\t{region_end}\tregion\n")
        tmp_path = tmp.name

    r = subprocess.run(['CrossMap', 'bed', chain_file, tmp_path, region_bed],
                       capture_output=True, text=True)
    os.unlink(tmp_path)

    region_start38, region_end38 = region_start, region_end
    if os.path.exists(region_bed):
        with open(region_bed) as f:
            line = f.readline().strip().split('\t')
            if len(line) >= 3:
                region_start38 = int(line[1]) + 1
                region_end38   = int(line[2])

    extract_path = os.path.join(output_folder, 'extract_positions.txt')

    plink_cmd = (
        f'/private/home/rsmerigl/plink2 '
        f'--vcf {vcf_file} '
        f'--chr {chr_val} '
        f'--from-bp {region_start38} '
        f'--to-bp {region_end38} '
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
    print(f"Submitted job {job_id}: {hit_label} (hg38 {region_start38}-{region_end38})")
