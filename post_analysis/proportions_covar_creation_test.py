"""
TEST: build the PROPORTIONS covariate file from the LAI output (GNomix or RFMix).

For each sample it computes the global ancestry proportions (fraction of the
window-haplotype assignments going to each ancestry, across all chromosomes),
then takes the base demographic covariates from input.covar (age, sex, BMI —
DROPPING the Global_PC* columns) and appends the proportions.

This replaces the PC-based covar with the ancestry-proportions covar that the
new pipeline uses (admixture mapping + conditional fine mapping condition on the
global ancestry proportions, with the reference ancestry dropped downstream).

Works for both GNomix and RFMix (both write MSP files); reads result_folder/
msp_folder/chr*.msp.

Output: result_folder/covar_file/covar_proportions.phe
        columns: IID age sex BMI <ancestries, alphabetical>
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import snputils as su

if len(sys.argv) != 3:
    print("Usage: python proportions_covar_creation_test.py <result_folder> <data_folder>")
    sys.exit(1)

result_folder = sys.argv[1]
data_folder   = sys.argv[2]

msp_glob   = os.path.join(result_folder, 'msp_folder', 'chr*.msp')
base_covar = os.path.join(data_folder, 'input.covar')
out_dir    = os.path.join(result_folder, 'covar_file')
out_covar  = os.path.join(out_dir, 'covar_proportions.phe')
os.makedirs(out_dir, exist_ok=True)


def to_code2name(amap):
    """Normalize an ancestry map {int->str} or {str->int} to {code(int): name}."""
    out = {}
    for k, v in amap.items():
        ks, vs = str(k).strip(), str(v).strip()
        if ks.isdigit() and not vs.isdigit():
            out[int(ks)] = vs
        elif vs.isdigit() and not ks.isdigit():
            out[int(vs)] = ks
    return out


# ── accumulate window-haplotype ancestry counts per sample across chromosomes ──
samples = None
code2name = None
counts = None          # DataFrame index=sample, columns=ancestry name

msp_files = sorted(glob.glob(msp_glob))
if not msp_files:
    raise FileNotFoundError(f'No MSP files found: {msp_glob}')

for msp in msp_files:
    lai = su.MSPReader(msp).read()
    if code2name is None:
        code2name = to_code2name(getattr(lai, 'ancestry_map'))
        samples = [str(s) for s in lai.samples]
        anc_names = sorted(code2name.values())
        counts = pd.DataFrame(0, index=samples, columns=anc_names, dtype=np.int64)

    arr = np.asarray(lai.lai, dtype=int)              # (n_windows, 2*n_samples)
    for cd, name in code2name.items():
        col_counts = (arr == cd).sum(axis=0)          # per haplotype column
        # haplotype columns are ordered [s0h0, s0h1, s1h0, s1h1, ...]
        samp_counts = col_counts.reshape(-1, 2).sum(axis=1)
        counts[name] += samp_counts

# ── proportions ────────────────────────────────────────────────────────────────
totals = counts.sum(axis=1)
prop = counts.div(totals, axis=0)                     # fraction per ancestry
prop = prop.reset_index().rename(columns={'index': 'IID'})
print(f'Ancestries found: {list(counts.columns)}')
print(f'Windows-haplotypes per sample (should be constant): '
      f'{totals.min()}–{totals.max()}')

# ── merge with base demographics (drop the PCs) ────────────────────────────────
cov = pd.read_csv(base_covar, sep='\t')
pc_cols = [c for c in cov.columns if c.startswith('Global_PC')]
cov = cov.drop(columns=pc_cols)
cov['IID'] = cov['IID'].astype(str)
prop['IID'] = prop['IID'].astype(str)

merged = cov.merge(prop, on='IID', how='inner')
# fixed column order: IID, age, sex, BMI, then ancestries (alphabetical)
demog = [c for c in ['IID', 'age', 'sex', 'BMI'] if c in merged.columns]
anc_cols = sorted(counts.columns)
merged = merged[demog + anc_cols]
merged.to_csv(out_covar, sep='\t', index=False)

print(f'\nSaved {len(merged)} samples -> {out_covar}')
print(f'Covar columns: {list(merged.columns)}')
# report the column indices so the downstream --covar-col-nums can be set
cols = list(merged.columns)
ref = 'EUR' if 'EUR' in anc_cols else anc_cols[-1]
ref_idx = cols.index(ref) + 1
keep_idx = [i + 1 for i, c in enumerate(cols) if c != 'IID' and c != ref]
print(f'Reference ancestry (dropped downstream): {ref} = column {ref_idx}')
print(f'Suggested --covar-col-nums (admixture): {",".join(str(i) for i in keep_idx)}  '
      f'(fine mapping: add LAI as column {len(cols)+1})')
