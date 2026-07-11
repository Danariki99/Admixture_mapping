"""
UKBB HLA-region ancestry-window counts.

Counts, for the classic MHC region (chr6:28,477,797-33,448,354, hg19), the number
of window x haplotype assignments per ancestry, computed BOTH on the complete
sample set (to validate against the old table) and on the KING-cutoff FILTERED
cohort.

Reads the chr6 RFMix .msp.tsv directly (no snputils) and streams it, so it only
fully parses the ~MHC windows.

Output: hla_region_ancestry_counts.csv  (rows: complete / filtered; cols: ancestries)
"""

import os
import numpy as np
import pandas as pd

MSP  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb/ukb_hap_chr6_v2_rfmix.msp.tsv'
KEEP = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
OUT  = os.path.join(os.path.dirname(__file__), 'hla_region_ancestry_counts.csv')

MHC_START, MHC_END = 25_726_063, 33_400_644          # extended MHC (xMHC), hg19
ANC = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']   # codes 0..7

# ── KING-cutoff keep IDs ──────────────────────────────────────────────────────
keep = set()
with open(KEEP) as f:
    for line in f:
        s = line.strip()
        if s and not s.startswith('#'):
            keep.add(s)
print(f'KING-cutoff keep IDs: {len(keep):,}')

# ── Stream the MSP ────────────────────────────────────────────────────────────
with open(MSP) as f:
    codes_line = f.readline()                        # subpopulation codes
    header = f.readline().rstrip('\n').split('\t')
    hap_cols = header[6:]                             # e.g. 2502845.0
    sample_ids = np.array([c.split('.')[0] for c in hap_cols])
    filt_idx = np.where(np.isin(sample_ids, list(keep)))[0]
    n_hap = len(hap_cols)
    print(f'Haplotype columns: {n_hap:,}  ({n_hap // 2:,} samples in MSP)')
    print(f'Filtered haplotype columns: {len(filt_idx):,}  ({len(filt_idx) // 2:,} samples)')

    full_counts = np.zeros(8, dtype=np.int64)
    filt_counts = np.zeros(8, dtype=np.int64)
    n_windows = 0

    for line in f:
        parts = line.split('\t', 6)
        spos, epos = int(parts[1]), int(parts[2])
        if not (spos < MHC_END and epos > MHC_START):   # overlap with MHC
            continue
        n_windows += 1
        codes = np.array(parts[6].split(), dtype=np.int8)
        if len(codes) != n_hap:
            raise ValueError(f'window {spos}-{epos}: {len(codes)} codes != {n_hap}')
        full_counts += np.bincount(codes, minlength=8)
        filt_counts += np.bincount(codes[filt_idx], minlength=8)

print(f'\nMHC windows counted: {n_windows}')

df = pd.DataFrame(
    [full_counts, filt_counts],
    index=['complete_samples', 'filtered_kingcutoff'],
    columns=ANC,
)
df['TOTAL'] = df[ANC].sum(axis=1)
df.to_csv(OUT)
print(f'\nSaved → {OUT}\n')

pd.set_option('display.width', 200)
print('=== HLA-region window x haplotype counts ===')
print(df.to_string())

print('\n=== proportions (%) per row ===')
prop = df[ANC].div(df['TOTAL'], axis=0) * 100
print(prop.round(3).to_string())
