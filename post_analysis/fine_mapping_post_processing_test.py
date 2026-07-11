import os
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from pyliftover import LiftOver

# TEST version of fine_mapping_post_processing.py — same conditional-analysis
# logic, parameterized on result_folder / data_folder.
#   BY (fdr_by) correction on the ADD and LAI p-values; a SNP is a candidate if
#   ADD is significant, LAI is NOT (the SNP explains the ancestry signal) and its
#   OR direction is concordant with the admixture-mapping window OR. Finally the
#   concordant candidates are lifted over hg19->hg38 into the Borzoi input VCF.

if len(sys.argv) < 3:
    print("Usage: python fine_mapping_post_processing_test.py <result_folder> <data_folder> [n_extensions 0-6]")
    sys.exit(1)

result_folder = sys.argv[1]
data_folder   = sys.argv[2]
N_EXT = int(sys.argv[3]) if len(sys.argv) > 3 else 0

fine_mapping_folder       = os.path.join(result_folder, 'fine_mapping_new')
fine_mapping_6wind_folder = os.path.join(result_folder, 'fine_mapping_new_6wind')
output_folder             = os.path.join(result_folder, 'fine_mapping_conditional_results')
output_file               = os.path.join(output_folder, 'fine_mapping_summary.tsv')
snp_list_file             = os.path.join(output_folder, 'fine_mapping_all_candidates.tsv')
# admixture-mapping GLM outputs (for the concordance direction)
admix_dir                 = os.path.join(result_folder, 'output')
BORZOI_VCF                = os.path.join(result_folder, 'borzoi_results', 'candidates_hg38_borzoi.vcf')

os.makedirs(output_folder, exist_ok=True)
ALPHA = 0.05


_admix_or_cache = {}
def admixture_window_or(ancestry, pheno, window_start):
    key = (ancestry, pheno, window_start)
    if key in _admix_or_cache:
        return _admix_or_cache[key]
    glm_file = os.path.join(admix_dir, f'output_ancestry_{ancestry}', pheno,
                            f'output.{pheno}.glm.logistic.hybrid')
    val = np.nan
    try:
        glm = pd.read_csv(glm_file, sep='\t')
        row = glm[(glm['TEST'] == 'ADD') & (glm['POS'] == window_start)]
        if not row.empty:
            val = pd.to_numeric(row['OR'].iloc[0], errors='coerce')
    except Exception:
        val = np.nan
    _admix_or_cache[key] = val
    return val


def load_hit_data(hit_path, allowed_windows=None):
    add_rows, lai_rows = [], []
    for filename in sorted(os.listdir(hit_path)):
        if not filename.endswith('.glm.logistic.hybrid'):
            continue
        base = filename.split('.')[0]
        file_parts = base.split('_')
        window_start, window_end = int(file_parts[-2]), int(file_parts[-1])
        if allowed_windows is not None and (window_start, window_end) not in allowed_windows:
            continue
        try:
            df = pd.read_csv(os.path.join(hit_path, filename), sep='\t')
            df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'PROV_REF', 'A1', 'OMITTED',
                          'A1_FREQ', 'FIRTH', 'TEST', 'OBS_CT', 'OR', 'LOG_OR_SE',
                          'L95', 'U95', 'Z_STAT', 'P', 'ERRCODE']
        except Exception as e:
            print(f'  Error reading {filename}: {e}')
            continue
        for col in ['P', 'OR', 'LOG_OR_SE', 'L95', 'U95']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df = df.dropna(subset=['P', 'OR'])
        df['window_start'], df['window_end'] = window_start, window_end
        add_rows.append(df[df['TEST'] == 'ADD'].copy())
        lai_rows.append(df[df['TEST'] == 'LAI'].copy())
    if not add_rows or not lai_rows:
        return None, None
    return pd.concat(add_rows, ignore_index=True), pd.concat(lai_rows, ignore_index=True)


def run_analysis(all_add, all_lai, label, ancestry, pheno):
    all_add = all_add.copy()
    all_add['P_BY'] = multipletests(all_add['P'].values, alpha=ALPHA, method='fdr_by')[1]
    all_lai = all_lai.copy()
    all_lai['P_BY'] = multipletests(all_lai['P'].values, alpha=ALPHA, method='fdr_by')[1]

    lai_cols = all_lai[['POS', 'ID', 'window_start', 'P', 'P_BY', 'OR', 'L95', 'U95']].rename(
        columns={'P': 'LAI_P', 'P_BY': 'LAI_P_BY', 'OR': 'LAI_OR', 'L95': 'LAI_L95', 'U95': 'LAI_U95'})
    merged = pd.merge(all_add, lai_cols, on=['POS', 'ID', 'window_start'], how='inner')

    merged['ADD_sig'] = merged['P_BY'] <= ALPHA
    merged['LAI_sig'] = merged['LAI_P_BY'] <= ALPHA
    candidates = merged[merged['ADD_sig'] & ~merged['LAI_sig']].copy()
    n_pre = len(candidates)
    if n_pre > 0:
        candidates['admixture_OR'] = candidates['window_start'].apply(
            lambda ws: admixture_window_or(ancestry, pheno, ws))
        candidates['concordant_direction'] = (
            np.sign(np.log(candidates['OR'])) == np.sign(np.log(candidates['admixture_OR'])))
        candidates = candidates[candidates['concordant_direction']].copy()

    stats = {'label': label, 'n_snps': len(merged), 'alpha': ALPHA,
             'n_add_sig': int(merged['ADD_sig'].sum()), 'n_add_not_sig': int((~merged['ADD_sig']).sum()),
             'n_lai_sig': int(merged['LAI_sig'].sum()), 'n_lai_not_sig': int((~merged['LAI_sig']).sum()),
             'n_pre_concordance': n_pre, 'n_discordant': n_pre - len(candidates),
             'n_candidates': len(candidates)}
    return merged, candidates, stats


all_candidates, summary_rows = [], []

if not os.path.isdir(fine_mapping_folder):
    print(f'No fine mapping folder: {fine_mapping_folder}')
    sys.exit(0)

for hit_dir in sorted(os.listdir(fine_mapping_folder)):
    hit_path = os.path.join(fine_mapping_folder, hit_dir)
    if not os.path.isdir(hit_path):
        continue
    parts = hit_dir.split('_')
    hit_ancestry, pheno, chr_num = parts[0], parts[1], parts[2].replace('chr', '')
    print(f'\n{hit_dir}')

    all_add, all_lai = load_hit_data(hit_path)
    if all_add is None:
        print('  Skipping: no data found')
        continue
    merged, candidates, stats = run_analysis(all_add, all_lai, 'sig windows only', hit_ancestry, pheno)

    if len(candidates) == 0 and N_EXT > 0 and os.path.isdir(os.path.join(fine_mapping_6wind_folder, hit_dir)):
        hit_6wind_path = os.path.join(fine_mapping_6wind_folder, hit_dir)
        ext_sorted = sorted((int(fn.split('.')[0].split('_')[-2]), int(fn.split('.')[0].split('_')[-1]))
                            for fn in os.listdir(hit_6wind_path) if fn.endswith('.glm.logistic.hybrid'))
        sig_windows = {(int(fn.split('.')[0].split('_')[-2]), int(fn.split('.')[0].split('_')[-1]))
                       for fn in os.listdir(hit_path) if fn.endswith('.glm.logistic.hybrid')}
        min_sig, max_sig = min(w[0] for w in sig_windows), max(w[1] for w in sig_windows)
        upstream   = [w for w in ext_sorted if w[1] <= min_sig][-N_EXT:]
        downstream = [w for w in ext_sorted if w[0] >= max_sig][:N_EXT]
        selected_ext = set(upstream) | set(downstream)
        all_add_ext, all_lai_ext = load_hit_data(hit_6wind_path, allowed_windows=selected_ext)
        if all_add_ext is not None:
            merged, candidates, stats = run_analysis(
                pd.concat([all_add, all_add_ext], ignore_index=True),
                pd.concat([all_lai, all_lai_ext], ignore_index=True),
                f'sig ± {N_EXT} windows', hit_ancestry, pheno)

    base_summary = {'hit': hit_dir, 'ancestry': hit_ancestry, 'pheno': pheno, 'chr': chr_num,
                    'analysis': stats['label'], 'alpha': ALPHA, 'n_snps': stats['n_snps'],
                    'n_add_sig': stats['n_add_sig'], 'n_add_not_sig': stats['n_add_not_sig'],
                    'n_lai_sig': stats['n_lai_sig'], 'n_lai_not_sig': stats['n_lai_not_sig'],
                    'n_pre_concordance': stats['n_pre_concordance'], 'n_discordant': stats['n_discordant'],
                    'n_candidates': len(candidates)}

    if len(candidates) == 0:
        print('  -> No candidates found')
        summary_rows.append(base_summary)
        continue

    candidates = candidates.copy()
    candidates['beta'] = np.log(candidates['OR'])
    candidates['hit'] = hit_dir
    candidates = candidates.sort_values('P_BY')
    hit_out = os.path.join(output_folder, hit_dir)
    os.makedirs(hit_out, exist_ok=True)
    out_cols = ['hit', 'ID', 'CHROM', 'POS', 'REF', 'ALT', 'A1', 'OR', 'beta', 'LOG_OR_SE',
                'L95', 'U95', 'P', 'P_BY', 'LAI_P', 'LAI_P_BY', 'LAI_OR', 'LAI_L95', 'LAI_U95',
                'admixture_OR', 'concordant_direction', 'window_start', 'window_end', 'OBS_CT']
    candidates[out_cols].to_csv(os.path.join(hit_out, f'{hit_dir}_candidates.tsv'), sep='\t', index=False)
    all_candidates.append(candidates[out_cols])

    top = candidates.loc[candidates['beta'].abs().idxmax()]
    base_summary.update({'top_snp': top['ID'], 'top_pos': top['POS'], 'top_OR': round(top['OR'], 4),
                         'top_beta': round(top['beta'], 4), 'top_P_BY': round(top['P_BY'], 6),
                         'top_LAI_P_BY': round(top['LAI_P_BY'], 6)})
    summary_rows.append(base_summary)
    print(f'  -> {len(candidates)} candidates | top: {top["ID"]} OR={top["OR"]:.3f}')

if summary_rows:
    pd.DataFrame(summary_rows).to_csv(output_file, sep='\t', index=False)
    print(f'\nSummary -> {output_file}')

if all_candidates:
    global_df = pd.concat(all_candidates, ignore_index=True)
    global_df.to_csv(snp_list_file, sep='\t', index=False)
    print(f'Candidates ({len(global_df)}) -> {snp_list_file}')

    # liftover hg19 -> hg38 for Borzoi
    os.makedirs(os.path.dirname(BORZOI_VCF), exist_ok=True)
    lo = LiftOver('hg19', 'hg38')
    unique_snps = global_df.drop_duplicates(subset='ID').copy()

    def liftover_pos(chrom, pos):
        r = lo.convert_coordinate(f'chr{chrom}', pos - 1)
        return (r[0][0].replace('chr', ''), r[0][1] + 1) if r else (None, None)

    unique_snps[['CHROM_hg38', 'POS_hg38']] = unique_snps.apply(
        lambda r: pd.Series(liftover_pos(r['CHROM'], r['POS'])), axis=1)
    vcf_snps = unique_snps.dropna(subset=['CHROM_hg38', 'POS_hg38']).copy()
    vcf_snps['CHROM_hg38'] = vcf_snps['CHROM_hg38'].apply(lambda x: str(int(float(x))))
    vcf_snps['POS_hg38'] = vcf_snps['POS_hg38'].astype(int)
    with open(BORZOI_VCF, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        for _, row in vcf_snps.iterrows():
            f.write(f"chr{row['CHROM_hg38']}\t{row['POS_hg38']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t.\t.\n")
    print(f'  {len(vcf_snps)} SNPs (hg38) -> {BORZOI_VCF}')
else:
    print('No candidates found across all hits.')
