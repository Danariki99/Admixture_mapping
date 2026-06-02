import os
import pandas as pd
from statsmodels.stats.multitest import multipletests

CANDIDATES_TSV = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
SUMMARY_TSV    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_summary.tsv'
PHENO_TABLE    = os.path.join(os.path.dirname(__file__), 'ukbb_v1.xlsx')
PHENO_SHEET    = 'first_batch'
OUT_XLSX       = os.path.join(os.path.dirname(__file__), 'table_fine_mapping.xlsx')


def load_pheno_map(path, sheet):
    df = pd.read_excel(path, sheet_name=sheet)
    return dict(zip(df['ID'], df['ID2']))


def build_table():
    df = pd.read_csv(CANDIDATES_TSV, sep='\t')
    summary = pd.read_csv(SUMMARY_TSV, sep='\t')

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

    # Bonferroni correction across all candidates
    _, _, _, p_bonf = multipletests(df['P'].values, alpha=0.05, method='bonferroni')
    df['P_bonferroni'] = p_bonf

    # OR with 95% CI as formatted strings
    df['OR (95% CI)']     = df.apply(lambda r: f"{r['OR']:.3f} ({r['L95']:.3f}–{r['U95']:.3f})", axis=1)
    df['LAI_OR (95% CI)'] = df.apply(lambda r: f"{r['LAI_OR']:.3f} ({r['LAI_L95']:.3f}–{r['LAI_U95']:.3f})"
                                     if pd.notna(r.get('LAI_OR')) else '', axis=1)

    # Concordance flag: ADD and LAI OR on same side of 1
    df['ADD_LAI_concordant'] = df.apply(
        lambda r: 'yes' if pd.notna(r.get('LAI_OR')) and ((r['OR'] > 1) == (r['LAI_OR'] > 1)) else 'no',
        axis=1
    )

    # Final column order — no hit-level constants (OBS_CT, n_snps, n_candidates)
    out = df[[
        'hit',
        'phenotype',
        'ancestry',
        'hit_chr',
        'ID',
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'A1',
        'OR (95% CI)',
        'OR',
        'beta',
        'LOG_OR_SE',
        'P',
        'P_BY',
        'P_bonferroni',
        'LAI_OR (95% CI)',
        'LAI_OR',
        'LAI_P',
        'LAI_P_BY',
        'ADD_LAI_concordant',
        'window_start',
        'window_end',
    ]].copy()

    out = out.rename(columns={
        'hit':               'Hit',
        'phenotype':         'Phenotype',
        'ancestry':          'Ancestry',
        'hit_chr':           'Chr',
        'ID':                'SNP ID',
        'CHROM':             'CHROM',
        'POS':               'POS (hg19)',
        'REF':               'REF',
        'ALT':               'ALT',
        'A1':                'Effect allele (A1)',
        'OR (95% CI)':       'ADD OR (95% CI)',
        'OR':                'ADD OR',
        'beta':              'ADD Beta',
        'LOG_OR_SE':         'ADD SE',
        'P':                 'ADD P (raw)',
        'P_BY':              'ADD P (BY-FDR)',
        'P_bonferroni':      'ADD P (Bonferroni)',
        'LAI_OR (95% CI)':   'LAI OR (95% CI)',
        'LAI_OR':            'LAI OR',
        'LAI_P':             'LAI P (raw)',
        'LAI_P_BY':          'LAI P (BY-FDR)',
        'ADD_LAI_concordant':'ADD/LAI concordant',
        'window_start':      'Window start',
        'window_end':        'Window end',
    })

    out = out.sort_values(['Hit', 'ADD P (raw)'])
    return out


def main():
    print('Building fine-mapping table...')
    out_df = build_table()
    out_df.to_excel(OUT_XLSX, index=False)
    print(f'  {len(out_df)} SNPs → {OUT_XLSX}')


if __name__ == '__main__':
    main()
