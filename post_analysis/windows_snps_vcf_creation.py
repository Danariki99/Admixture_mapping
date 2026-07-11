"""
Builds the hg38 VCF of ALL SNPs tested in the significant fine-mapping windows
(not just the candidates) — for running Borzoi on the full window content.

Reads every ADD row from fine_mapping_new/*/*.glm.logistic.hybrid, deduplicates
by SNP ID, lifts over hg19 -> hg38, and writes windows_all_hg38.vcf.

Output: borzoi_results/windows_all_hg38.vcf
"""

import os
import glob
import subprocess
import pandas as pd
from pyliftover import LiftOver

FINE_MAP_DIR = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
OUT_VCF      = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/windows_all_hg38.vcf'

# ── 1. Unique SNPs (ID -> chrom, pos, ref, alt) from all ADD rows ──────────────
print('Extracting unique ADD SNPs from fine_mapping_new ...')
files = glob.glob(os.path.join(FINE_MAP_DIR, '*', '*.glm.logistic.hybrid'))
# awk: TEST is col 11, print CHROM POS ID REF ALT for ADD rows; unique by ID handled below
awk = r"""awk -F'\t' '$11=="ADD"{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'"""
res = subprocess.run(f'{awk} {" ".join(files)}', shell=True, capture_output=True, text=True)
rows = [l.split('\t') for l in res.stdout.strip().split('\n') if l]
df = pd.DataFrame(rows, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT'])
df = df.drop_duplicates(subset='ID').reset_index(drop=True)
df['POS'] = df['POS'].astype(int)
print(f'  Unique SNPs: {len(df):,}')

# ── 2. Liftover hg19 -> hg38 ───────────────────────────────────────────────────
print('Lifting over hg19 -> hg38 ...')
lo = LiftOver('hg19', 'hg38')

def liftover_pos(chrom, pos):
    r = lo.convert_coordinate(f'chr{chrom}', pos - 1)   # pyliftover is 0-based
    if r:
        return r[0][0].replace('chr', ''), r[0][1] + 1
    return None, None

df[['CHROM_hg38', 'POS_hg38']] = df.apply(
    lambda r: pd.Series(liftover_pos(r['CHROM'], r['POS'])), axis=1)

n_fail = df['POS_hg38'].isna().sum()
if n_fail:
    print(f'  WARNING: {n_fail} SNPs failed liftover, skipped')
df = df.dropna(subset=['CHROM_hg38', 'POS_hg38']).copy()
df['CHROM_hg38'] = df['CHROM_hg38'].apply(lambda x: str(int(float(x))))
df['POS_hg38'] = df['POS_hg38'].astype(int)

# ── 3. Write VCF ───────────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUT_VCF), exist_ok=True)
with open(OUT_VCF, 'w') as f:
    f.write('##fileformat=VCFv4.2\n')
    for _, r in df.iterrows():
        f.write(f"chr{r['CHROM_hg38']}\t{r['POS_hg38']}\t{r['ID']}\t{r['REF']}\t{r['ALT']}\t.\t.\n")

print(f'\n{len(df):,} SNPs -> {OUT_VCF}')
