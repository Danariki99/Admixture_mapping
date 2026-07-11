"""
Ancestry-proportion comparison table (KING-cutoff FILTERED cohort):
  - Global      : genome-wide window x haplotype ancestry counts
  - HLA (xMHC)  : counts inside the extended MHC (from hla_region_ancestry_counts.csv)
  - Non-HLA     : Global - HLA

Global counts come from the per-sample per-chromosome final_counts_chr*.csv,
summed over the KING-cutoff samples and over all autosomes. HLA counts are the
'filtered_kingcutoff' row of hla_region_ancestry_counts.csv (extended MHC, hg19).

Outputs:
  hla_global_comparison_counts.csv        (Global / HLA / Non-HLA, raw counts)
  hla_global_comparison_proportions.csv   (same, row-normalised %)
"""

import os
import numpy as np
import pandas as pd
import openpyxl

COUNTS_DIR = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files_old/ukbb/counts'
KEEP       = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
HLA_CSV    = os.path.join(os.path.dirname(__file__), 'hla_region_ancestry_counts.csv')
OUT_COUNTS = os.path.join(os.path.dirname(__file__), 'hla_global_comparison_counts.csv')
OUT_PROP   = os.path.join(os.path.dirname(__file__), 'hla_global_comparison_proportions.csv')

ANC = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']

# ── KING-cutoff keep IDs ──────────────────────────────────────────────────────
keep = pd.read_csv(KEEP, comment=None)
keep = pd.read_csv(KEEP); keep.columns = ['IID']
keep_ids = set(keep['IID'].astype(str))
print(f'KING-cutoff keep IDs: {len(keep_ids):,}')

# ── Genome-wide global counts (filtered), summing all autosomes ───────────────
global_counts = pd.Series(0, index=ANC, dtype=np.int64)
matched = None
for chrom in range(1, 23):
    f = os.path.join(COUNTS_DIR, f'final_counts_chr{chrom}.csv')
    df = pd.read_csv(f)
    df['#IID'] = df['#IID'].astype(str)
    df = df[df['#IID'].isin(keep_ids)]
    if matched is None:
        matched = len(df)
        print(f'  chr1 matched samples: {matched:,}')
    global_counts += df[ANC].sum()
print(f'Global counts computed over {matched:,} filtered samples\n')

# ── HLA (extended MHC) counts, filtered ───────────────────────────────────────
hla = pd.read_csv(HLA_CSV, index_col=0)
hla_counts = hla.loc['filtered_kingcutoff', ANC].astype(np.int64)

# ── Non-HLA = Global - HLA ────────────────────────────────────────────────────
nonhla_counts = global_counts - hla_counts

table = pd.DataFrame(
    [global_counts, hla_counts, nonhla_counts],
    index=['Global', 'HLA_xMHC', 'Non_HLA'],
)[ANC]
table['TOTAL'] = table.sum(axis=1)
table.to_csv(OUT_COUNTS)
table.to_excel(OUT_COUNTS.replace('.csv', '.xlsx'), index=True, index_label='Region') 

prop = table[ANC].div(table['TOTAL'], axis=0) * 100
prop.to_csv(OUT_PROP)
prop.to_excel(OUT_PROP.replace('.csv', '.xlsx'), index=True, index_label='Region')

pd.set_option('display.width', 220)
print('=== counts ===')
print(table.to_string())
print('\n=== proportions (%) ===')
print(prop.round(3).to_string())
print(f'\nSaved → {OUT_COUNTS}')
print(f'Saved → {OUT_PROP}')

