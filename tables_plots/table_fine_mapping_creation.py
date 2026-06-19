import os
import pandas as pd

CANDIDATES_TSV = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
CAUSAL_CANDIDATES_TSV = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_causal_candidates.tsv'
ANNOTATED_TSV  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_candidates_annotated.tsv'
SUMMARY_TSV    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_summary.tsv'
PHENO_TABLE    = os.path.join(os.path.dirname(__file__), 'ukbb_v1.xlsx')
PHENO_SHEET    = 'first_batch'
OUT_XLSX       = os.path.join(os.path.dirname(__file__), 'table_fine_mapping.xlsx')
OUT_CAUSAL_XLSX = os.path.join(os.path.dirname(__file__), 'table_fine_mapping_causal.xlsx')


def load_pheno_map(path, sheet):
    df = pd.read_excel(path, sheet_name=sheet)
    return dict(zip(df['ID'], df['ID2']))


# ── Causal-variant table ──────────────────────────────────────────────────────
# The causal SNP of each locus (the variant with the largest effect for that
# phenotype) is taken from CAUSAL_CANDIDATES_TSV. Several loci share the same
# causal SNP, so the set collapses to a handful of unique variants; the phenotype
# hit(s) a SNP is causal for are then listed together on its row.

# Per-ancestry alternate-allele frequencies (gnomAD) for the causal variants, as
# decimal fractions. Only the frequency of the involved ancestry is reported in
# the table, formatted as a percentage.
CAUSAL_FREQ = {
    'rs61737338': {'AFR': 0.06781, 'EAS': 0.00002230, 'EUR': 0.005570, 'AMR': 0.02233, 'SAS': 0.0006041, 'WAS': 0.01255},
    'rs3177928':  {'AFR': 0.08487, 'EAS': 0.04213,    'EUR': 0.1484,   'AMR': 0.1075,  'SAS': 0.1405,    'WAS': 0.05782},
    'rs6906021':  {'AFR': 0.4967,  'EAS': 0.5089,     'EUR': 0.4622,   'AMR': 0.4602,  'SAS': 0.5692,    'WAS': 0.6575},
    'rs9274569':  {'AFR': 0.3948,  'EAS': 0.2868,     'EUR': 0.4098,   'AMR': 0.3540,  'SAS': 0.3614,    'WAS': 0.5476},
    'rs28732226': {'SAS': 0.1385},
}

FREQ_COL = 'Allele frequency (involved ancestry, %)'


def _load_annotated():
    """Load all fine-mapping candidates, merge VEP annotation, parse the hit
    name, and add formatted OR/CI strings plus the ADD/LAI concordance flag.
    Returns one row per (hit, SNP) candidate."""
    df = pd.read_csv(CANDIDATES_TSV, sep='\t')

    # Merge VEP annotation (unique per SNP ID)
    if os.path.exists(ANNOTATED_TSV):
        ann = pd.read_csv(ANNOTATED_TSV, sep='\t')[['ID', 'category', 'gene_symbol']].drop_duplicates('ID')
        ann['gene_symbol'] = ann['gene_symbol'].fillna('None')
        df = df.merge(ann, on='ID', how='left')
    else:
        df['category']    = 'N/A'
        df['gene_symbol'] = 'N/A'

    # Phenotype name map
    pheno_map = load_pheno_map(PHENO_TABLE, PHENO_SHEET) if os.path.exists(PHENO_TABLE) else {}

    # Parse ancestry, phenotype ID, chromosome from hit name (e.g. AFR_HC219_chr6)
    def parse_hit(hit):
        parts = hit.split('_')
        ancestry = parts[0]
        chrom    = parts[-1]
        pheno_id = '_'.join(parts[1:-1])
        return ancestry, pheno_id, chrom

    df[['ancestry', 'pheno_id', 'hit_chr']] = df['hit'].apply(
        lambda h: pd.Series(parse_hit(h))
    )
    df['phenotype'] = df['pheno_id'].map(pheno_map).fillna(df['pheno_id'])
    df['phenotype'] = df['phenotype'].apply(lambda s: ' '.join(s.split('_')[1:]) if '_' in s else s)

    # OR with 95% CI as formatted strings (2 decimal places)
    df['OR (95% CI)']     = df.apply(lambda r: f"{r['OR']:.2f} ({r['L95']:.2f}–{r['U95']:.2f})", axis=1)
    df['LAI_OR (95% CI)'] = df.apply(lambda r: f"{r['LAI_OR']:.2f} ({r['LAI_L95']:.2f}–{r['LAI_U95']:.2f})"
                                     if pd.notna(r.get('LAI_OR')) else '', axis=1)

    # Concordance flag: ADD and LAI OR on same side of 1
    df['ADD_LAI_concordant'] = df.apply(
        lambda r: 'yes' if pd.notna(r.get('LAI_OR')) and ((r['OR'] > 1) == (r['LAI_OR'] > 1)) else 'no',
        axis=1
    )
    return df


RENAME_MAP = {
    'phenotype':         'Phenotype',
    'ancestry':          'Ancestry',
    'hit_chr':           'Chr',
    'ID':                'SNP ID',
    'CHROM':             'CHROM',
    'POS':               'POS (hg19)',
    'REF':               'REF',
    'ALT':               'ALT',
    'A1':                'Effect allele (A1)',
    'gene_symbol':       'Gene',
    'category':          'Variant type',
    'OR (95% CI)':       'ADD OR (95% CI)',
    'P':                 'ADD P (raw)',
    'P_BY':              'ADD P (BY-FDR)',
    'LAI_OR (95% CI)':   'LAI OR (95% CI)',
    'LAI_P':             'LAI P (raw)',
    'LAI_P_BY':          'LAI P (BY-FDR)',
    'ADD_LAI_concordant':'ADD/LAI concordant',
    'window_start':      'Window start',
    'window_end':        'Window end',
}


def build_table():
    df = _load_annotated()

    # Final column order
    out = df[[
        'phenotype',
        'ancestry',
        'hit_chr',
        'ID',
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'A1',
        'gene_symbol',
        'category',
        'OR (95% CI)',
        'P',
        'P_BY',
        'LAI_OR (95% CI)',
        'LAI_P',
        'LAI_P_BY',
        'ADD_LAI_concordant',
        'window_start',
        'window_end',
    ]].copy()

    out = out.rename(columns=RENAME_MAP)

    out = out.sort_values(['Phenotype', 'Ancestry', 'Chr', 'ADD P (raw)'])
    return out


def build_causal_table():
    """Reduce the candidate table to the causal variants reported in the paper:
    one row per unique causal SNP, with every phenotype it is causal for listed
    alongside, plus the allele frequency of the involved ancestry.

    The causal SNP of each phenotype hit is read from CAUSAL_CANDIDATES_TSV
    (largest effect per locus); rows are then collapsed to unique SNPs so a SNP
    that is causal for several phenotypes appears once with the phenotypes joined."""
    df = _load_annotated()

    # Causal SNP per hit (one row per hit)
    causal = pd.read_csv(CAUSAL_CANDIDATES_TSV, sep='\t')[['hit', 'candidate_snp']].drop_duplicates()

    # Pull the full annotated stats for each (hit, causal SNP)
    sel = df.merge(causal, left_on=['hit', 'ID'],
                   right_on=['hit', 'candidate_snp'], how='inner')

    missing = set(causal['candidate_snp']) - set(sel['ID'])
    if missing:
        print(f"  WARNING: causal SNP(s) not found in candidates: {sorted(missing)}")

    rows = []
    for snp, g in sel.groupby('ID'):
        ancestry = g['ancestry'].iloc[0]
        # Representative row = the hit where the SNP is most significant
        rep = g.loc[g['P'].idxmin()].copy()
        # Phenotype(s) the SNP is causal for, ordered by significance (deduped)
        phenos = list(dict.fromkeys(g.sort_values('P')['phenotype']))
        rep['phenotype'] = ', '.join(phenos)
        rep['ancestry']  = ancestry
        freq             = CAUSAL_FREQ.get(snp, {}).get(ancestry)
        rep[FREQ_COL]    = f"{freq * 100:.2f}%" if freq is not None else None
        rows.append(rep)

    out = pd.DataFrame(rows)
    out = out[[
        'phenotype', 'ancestry', 'hit_chr', 'ID', 'CHROM', 'POS', 'REF', 'ALT', 'A1',
        'gene_symbol', 'category', FREQ_COL,
        'OR (95% CI)', 'P', 'P_BY',
        'LAI_OR (95% CI)', 'LAI_P', 'LAI_P_BY', 'ADD_LAI_concordant',
        'window_start', 'window_end',
    ]].rename(columns=RENAME_MAP)

    out = out.sort_values(['Ancestry', 'ADD P (raw)'])
    return out


def main():
    print('Building fine-mapping table...')
    out_df = build_table()
    out_df.to_excel(OUT_XLSX, index=False)
    print(f'  {len(out_df)} SNPs → {OUT_XLSX}')

    print('Building causal-variant table (one row per unique causal SNP)...')
    causal_df = build_causal_table()
    causal_df.to_excel(OUT_CAUSAL_XLSX, index=False)
    print(f'  {len(causal_df)} causal SNPs → {OUT_CAUSAL_XLSX}')


if __name__ == '__main__':
    main()
