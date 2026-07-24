import os
import sys
import glob
import subprocess
import pandas as pd
import pysam
from pyliftover import LiftOver

# TEST version of windows_snps_vcf_creation.py.
# Builds the hg38 VCF of ALL SNPs tested in the fine-mapping windows (not only the
# candidates) — the input for running Borzoi on the full window content.
# Parameterized on result_folder.

if len(sys.argv) < 2:
    print("Usage: python windows_snps_vcf_creation_test.py <result_folder>")
    sys.exit(1)

result_folder = sys.argv[1]
FINE_MAP_DIR = os.path.join(result_folder, 'fine_mapping_new')
OUT_VCF      = os.path.join(result_folder, 'borzoi_results', 'windows_all_hg38.vcf')

files = glob.glob(os.path.join(FINE_MAP_DIR, '*', '*.glm.logistic.hybrid'))
if not files:
    print(f'No fine-mapping glm files under {FINE_MAP_DIR}. Nothing to do.')
    sys.exit(0)

print('Extracting unique ADD SNPs from fine_mapping_new ...')
awk = r"""awk -F'\t' '$11=="ADD"{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'"""
res = subprocess.run(f'{awk} {" ".join(files)}', shell=True, capture_output=True, text=True)
rows = [l.split('\t') for l in res.stdout.strip().split('\n') if l]
df = pd.DataFrame(rows, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT']).drop_duplicates('ID').reset_index(drop=True)
df['POS'] = df['POS'].astype(int)
print(f'  Unique SNPs: {len(df):,}')

print('Lifting over hg19 -> hg38 ...')
lo = LiftOver('hg19', 'hg38')

def liftover_pos(chrom, pos):
    r = lo.convert_coordinate(f'chr{chrom}', pos - 1)
    return (r[0][0].replace('chr', ''), r[0][1] + 1) if r else (None, None)

df[['CHROM_hg38', 'POS_hg38']] = df.apply(lambda r: pd.Series(liftover_pos(r['CHROM'], r['POS'])), axis=1)
n_fail = df['POS_hg38'].isna().sum()
if n_fail:
    print(f'  WARNING: {n_fail} SNPs failed liftover, skipped')
df = df.dropna(subset=['CHROM_hg38', 'POS_hg38']).copy()
df['CHROM_hg38'] = df['CHROM_hg38'].apply(lambda x: str(int(float(x))))
df['POS_hg38'] = df['POS_hg38'].astype(int)

# Drop SNPs whose REF does not match the hg38 reference (liftover strand flips):
# Borzoi builds an EMPTY sequence for them and CRASHES the whole run, leaving the
# pre-allocated sed.h5 full of zeros.
HG38_FA = os.environ.get(
    'BORZOI_HG38_FA',
    '/private/home/rsmerigl/codes/cleaned_codes/borzoi/examples/hg38/assembly/ucsc/hg38.fa')
_fa = pysam.Fastafile(HG38_FA)

def _ref_matches_hg38(r):
    base = _fa.fetch(f"chr{r['CHROM_hg38']}", r['POS_hg38'] - 1,
                     r['POS_hg38'] - 1 + len(r['REF'])).upper()
    return base == r['REF'].upper()

_ok = df.apply(_ref_matches_hg38, axis=1)
if (~_ok).any():
    print(f'  Dropping {int((~_ok).sum())} SNPs whose REF does not match hg38 (strand flips)')
df = df[_ok].copy()

# sort by (chromosome, position): keeps chromosomes contiguous and lets Borzoi
# cluster neighbouring SNPs
df = df.sort_values(['CHROM_hg38', 'POS_hg38'], key=lambda c: c.astype(int) if c.name == 'CHROM_hg38' else c)

os.makedirs(os.path.dirname(OUT_VCF), exist_ok=True)
with open(OUT_VCF, 'w') as f:
    f.write('##fileformat=VCFv4.2\n')
    for _, r in df.iterrows():
        f.write(f"chr{r['CHROM_hg38']}\t{r['POS_hg38']}\t{r['ID']}\t{r['REF']}\t{r['ALT']}\t.\t.\n")
print(f'\n{len(df):,} SNPs -> {OUT_VCF}')
