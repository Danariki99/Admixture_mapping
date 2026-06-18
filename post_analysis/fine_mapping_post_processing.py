import os
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from pyliftover import LiftOver

if len(sys.argv) != 2 or not sys.argv[1].isdigit() or not (0 <= int(sys.argv[1]) <= 6):
    print("Usage: python fine_mapping_post_processing.py <n_extensions (0-6)>")
    print("  0 = significant windows only (no extension)")
    sys.exit(1)

N_EXT = int(sys.argv[1])

fine_mapping_folder       = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
fine_mapping_6wind_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new_6wind'
output_folder             = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results'
output_file               = os.path.join(output_folder, 'fine_mapping_summary.tsv')
snp_list_file             = os.path.join(output_folder, 'fine_mapping_all_candidates.tsv')
admix_dir                 = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb'

os.makedirs(output_folder, exist_ok=True)


ALPHA = 0.05


# Window-level admixture-mapping OR (ADD test) — for concordance filtering
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

        base       = filename.split('.')[0]
        file_parts = base.split('_')
        window_start = int(file_parts[-2])
        window_end   = int(file_parts[-1])

        if allowed_windows is not None and (window_start, window_end) not in allowed_windows:
            continue

        filepath = os.path.join(hit_path, filename)
        try:
            df = pd.read_csv(filepath, sep='\t')
            df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'PROV_REF', 'A1', 'OMITTED',
                          'A1_FREQ', 'FIRTH', 'TEST', 'OBS_CT', 'OR', 'LOG_OR_SE',
                          'L95', 'U95', 'Z_STAT', 'P', 'ERRCODE']
        except Exception as e:
            print(f'  Error reading {filename}: {e}')
            continue

        for col in ['P', 'OR', 'LOG_OR_SE', 'L95', 'U95']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df = df.dropna(subset=['P', 'OR'])
        df['window_start'] = window_start
        df['window_end']   = window_end

        add_rows.append(df[df['TEST'] == 'ADD'].copy())
        lai_rows.append(df[df['TEST'] == 'LAI'].copy())

    if not add_rows or not lai_rows:
        return None, None

    return pd.concat(add_rows, ignore_index=True), pd.concat(lai_rows, ignore_index=True)


def run_analysis(all_add, all_lai, label, ancestry, pheno):
    """
    BY on ADD and LAI p-values; fixed threshold ALPHA for significance.
    Candidates: ADD_P_BY <= ALPHA AND LAI_P_BY > ALPHA AND OR direction concordant
    with the admixture-mapping window-level OR.
    """
    _, add_by, _, _ = multipletests(all_add['P'].values, alpha=ALPHA, method='fdr_by')
    all_add = all_add.copy()
    all_add['P_BY'] = add_by

    _, lai_by, _, _ = multipletests(all_lai['P'].values, alpha=ALPHA, method='fdr_by')
    all_lai = all_lai.copy()
    all_lai['P_BY'] = lai_by

    lai_cols = all_lai[['POS', 'ID', 'window_start', 'P', 'P_BY', 'OR', 'L95', 'U95']].rename(
        columns={'P': 'LAI_P', 'P_BY': 'LAI_P_BY', 'OR': 'LAI_OR', 'L95': 'LAI_L95', 'U95': 'LAI_U95'}
    )
    merged = pd.merge(all_add, lai_cols, on=['POS', 'ID', 'window_start'], how='inner')

    merged['ADD_sig'] = merged['P_BY'] <= ALPHA
    merged['LAI_sig'] = merged['LAI_P_BY'] <= ALPHA

    candidates = merged[merged['ADD_sig'] & ~merged['LAI_sig']].copy()
    n_pre_concordance = len(candidates)

    # Keep only candidates whose OR direction is concordant with the
    # admixture-mapping window-level OR direction.
    if n_pre_concordance > 0:
        candidates['admixture_OR'] = candidates['window_start'].apply(
            lambda ws: admixture_window_or(ancestry, pheno, ws)
        )
        candidates['concordant_direction'] = (
            np.sign(np.log(candidates['OR'])) == np.sign(np.log(candidates['admixture_OR']))
        )
        candidates = candidates[candidates['concordant_direction']].copy()

    stats = {
        'label':              label,
        'n_snps':             len(merged),
        'alpha':              ALPHA,
        'n_add_sig':          int(merged['ADD_sig'].sum()),
        'n_add_not_sig':      int((~merged['ADD_sig']).sum()),
        'n_lai_sig':          int(merged['LAI_sig'].sum()),
        'n_lai_not_sig':      int((~merged['LAI_sig']).sum()),
        'n_pre_concordance':  n_pre_concordance,
        'n_discordant':       n_pre_concordance - len(candidates),
        'n_candidates':       len(candidates),
    }
    return merged, candidates, stats


def print_stats(hit_dir, stats):
    print(f"  [{stats['label']}]  SNPs: {stats['n_snps']}  (BY threshold: {stats['alpha']})")
    print(f"  ADD sig: {stats['n_add_sig']}  |  ADD not sig: {stats['n_add_not_sig']}")
    print(f"  LAI sig: {stats['n_lai_sig']}  |  LAI not sig: {stats['n_lai_not_sig']}")
    print(f"  Candidates (ADD sig + LAI not sig): {stats['n_pre_concordance']}  "
          f"-> concordant: {stats['n_candidates']}  (discordant removed: {stats['n_discordant']})")


all_candidates = []
summary_rows   = []

for hit_dir in sorted(os.listdir(fine_mapping_folder)):
    hit_path = os.path.join(fine_mapping_folder, hit_dir)
    if not os.path.isdir(hit_path):
        continue

    parts        = hit_dir.split('_')
    hit_ancestry = parts[0]
    pheno        = parts[1]
    chr_num      = parts[2].replace('chr', '')

    print(f'\n{hit_dir}')

    # Step 1: significant windows only
    all_add, all_lai = load_hit_data(hit_path)
    if all_add is None:
        print('  Skipping: no data found')
        continue

    merged, candidates, stats = run_analysis(all_add, all_lai, label='sig windows only',
                                             ancestry=hit_ancestry, pheno=pheno)
    print_stats(hit_dir, stats)

    # Step 2: extend ±N_EXT if no candidates
    if len(candidates) == 0 and N_EXT > 0:
        print(f'  No candidates in sig windows — trying ±{N_EXT} extension windows...')

        hit_6wind_path = os.path.join(fine_mapping_6wind_folder, hit_dir)
        if not os.path.isdir(hit_6wind_path):
            print(f'  Extension folder not found: {hit_6wind_path}')
        else:
            ext_windows_sorted = []
            for fn in sorted(os.listdir(hit_6wind_path)):
                if not fn.endswith('.glm.logistic.hybrid'):
                    continue
                fp = fn.split('.')[0].split('_')
                ext_windows_sorted.append((int(fp[-2]), int(fp[-1])))
            ext_windows_sorted.sort()

            sig_windows = set()
            for fn in os.listdir(hit_path):
                if not fn.endswith('.glm.logistic.hybrid'):
                    continue
                fp = fn.split('.')[0].split('_')
                sig_windows.add((int(fp[-2]), int(fp[-1])))

            min_sig = min(w[0] for w in sig_windows)
            max_sig = max(w[1] for w in sig_windows)

            upstream     = [w for w in ext_windows_sorted if w[1] <= min_sig][-N_EXT:] if N_EXT > 0 else []
            downstream   = [w for w in ext_windows_sorted if w[0] >= max_sig][:N_EXT]  if N_EXT > 0 else []
            selected_ext = set(upstream) | set(downstream)

            print(f'  Using {len(upstream)} upstream + {len(sig_windows)} sig + {len(downstream)} downstream windows')

            all_add_ext, all_lai_ext = load_hit_data(hit_6wind_path, allowed_windows=selected_ext)
            if all_add_ext is not None:
                all_add_ext = pd.concat([all_add, all_add_ext], ignore_index=True)
                all_lai_ext = pd.concat([all_lai, all_lai_ext], ignore_index=True)
                merged, candidates, stats = run_analysis(
                    all_add_ext, all_lai_ext, label=f'sig ± {N_EXT} windows',
                    ancestry=hit_ancestry, pheno=pheno
                )
                print_stats(hit_dir, stats)

    if len(candidates) == 0:
        print('  -> No candidates found')
        summary_rows.append({
            'hit':           hit_dir,
            'ancestry':      hit_ancestry,
            'pheno':         pheno,
            'chr':           chr_num,
            'analysis':      stats['label'],
            'alpha':         ALPHA,
            'n_snps':        stats['n_snps'],
            'n_add_sig':     stats['n_add_sig'],
            'n_add_not_sig': stats['n_add_not_sig'],
            'n_lai_sig':     stats['n_lai_sig'],
            'n_lai_not_sig': stats['n_lai_not_sig'],
            'n_pre_concordance': stats['n_pre_concordance'],
            'n_discordant':  stats['n_discordant'],
            'n_candidates':  0,
        })
        continue

    # Annotate and save per-hit candidate list
    candidates = candidates.copy()
    candidates['beta'] = np.log(candidates['OR'])
    candidates['hit']  = hit_dir
    candidates = candidates.sort_values('P_BY')

    hit_out = os.path.join(output_folder, hit_dir)
    os.makedirs(hit_out, exist_ok=True)
    cand_file = os.path.join(hit_out, f'{hit_dir}_candidates.tsv')
    out_cols = ['hit', 'ID', 'CHROM', 'POS', 'REF', 'ALT', 'A1',
                'OR', 'beta', 'LOG_OR_SE', 'L95', 'U95',
                'P', 'P_BY', 'LAI_P', 'LAI_P_BY', 'LAI_OR', 'LAI_L95', 'LAI_U95',
                'admixture_OR', 'concordant_direction',
                'window_start', 'window_end', 'OBS_CT']
    candidates[out_cols].to_csv(cand_file, sep='\t', index=False)

    all_candidates.append(candidates[out_cols])

    top = candidates.loc[candidates['beta'].abs().idxmax()]
    print(f'  -> {len(candidates)} candidates | top effect: {top["ID"]} chr{top["CHROM"]}:{top["POS"]} '
          f'OR={top["OR"]:.3f} beta={top["beta"]:.3f} P_BY={top["P_BY"]:.4g} LAI_P_BY={top["LAI_P_BY"]:.4g}')
    print(f'  Saved -> {cand_file}')

    summary_rows.append({
        'hit':              hit_dir,
        'ancestry':         hit_ancestry,
        'pheno':            pheno,
        'chr':              chr_num,
        'analysis':         stats['label'],
        'alpha':            ALPHA,
        'n_snps':           stats['n_snps'],
        'n_add_sig':        stats['n_add_sig'],
        'n_add_not_sig':    stats['n_add_not_sig'],
        'n_lai_sig':        stats['n_lai_sig'],
        'n_lai_not_sig':    stats['n_lai_not_sig'],
        'n_pre_concordance': stats['n_pre_concordance'],
        'n_discordant':     stats['n_discordant'],
        'n_candidates':     len(candidates),
        'top_snp':          top['ID'],
        'top_pos':          top['POS'],
        'top_OR':           round(top['OR'], 4),
        'top_beta':         round(top['beta'], 4),
        'top_P_BY':         round(top['P_BY'], 6),
        'top_LAI_P_BY':     round(top['LAI_P_BY'], 6),
    })

print(f'\n{"="*60}')

if summary_rows:
    pd.DataFrame(summary_rows).to_csv(output_file, sep='\t', index=False)
    print(f'Summary saved -> {output_file}')

if all_candidates:
    global_df = pd.concat(all_candidates, ignore_index=True)
    global_df.to_csv(snp_list_file, sep='\t', index=False)
    print(f'Global SNP list ({len(global_df)} candidates) -> {snp_list_file}')
else:
    print('No candidates found across all hits.')

# ── Liftover hg19 → hg38 and write VCF for Borzoi ────────────────────────────
BORZOI_VCF = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/candidates_hg38_borzoi.vcf'

if all_candidates:
    print(f'\nLifting over candidates hg19 → hg38 for Borzoi...')
    lo = LiftOver('hg19', 'hg38')

    # Unique SNPs only (same SNP can appear in multiple hits)
    unique_snps = global_df.drop_duplicates(subset='ID').copy()

    def liftover_pos(chrom, pos):
        result = lo.convert_coordinate(f'chr{chrom}', pos - 1)  # pyliftover is 0-based
        if result and len(result) > 0:
            chrom_hg38 = result[0][0].replace('chr', '')
            pos_hg38   = result[0][1] + 1  # back to 1-based
            return chrom_hg38, pos_hg38
        return None, None

    unique_snps[['CHROM_hg38', 'POS_hg38']] = unique_snps.apply(
        lambda r: pd.Series(liftover_pos(r['CHROM'], r['POS'])), axis=1
    )

    n_failed = unique_snps['POS_hg38'].isna().sum()
    if n_failed > 0:
        print(f'  WARNING: {n_failed} SNPs failed liftover and will be skipped')

    vcf_snps = unique_snps.dropna(subset=['CHROM_hg38', 'POS_hg38']).copy()
    vcf_snps['CHROM_hg38'] = vcf_snps['CHROM_hg38'].apply(lambda x: str(int(float(x))))
    vcf_snps['POS_hg38']   = vcf_snps['POS_hg38'].astype(int)

    with open(BORZOI_VCF, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        for _, row in vcf_snps.iterrows():
            f.write(f"chr{row['CHROM_hg38']}\t{row['POS_hg38']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t.\t.\n")

    print(f'  {len(vcf_snps)} SNPs → {BORZOI_VCF}')
    print(f'  Run: diff {BORZOI_VCF} /private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/candidates_hg38.vcf')
