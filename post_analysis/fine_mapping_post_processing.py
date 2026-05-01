import os
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

if len(sys.argv) != 2 or not sys.argv[1].isdigit() or not (0 <= int(sys.argv[1]) <= 6):
    print("Usage: python fine_mapping_post_processing.py <n_extensions (0-6)>")
    print("  0 = significant windows only (no extension)")
    sys.exit(1)

N_EXT = int(sys.argv[1])

fine_mapping_folder      = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
fine_mapping_6wind_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new_6wind'
wind_folder              = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
output_file              = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new/fine_mapping_results.tsv'


def load_hit_data(hit_path, allowed_windows=None):
    """
    Read all .glm.logistic.hybrid files in hit_path.
    If allowed_windows is a set of (start, end) tuples, only load those windows.
    Returns (all_add, all_lai) DataFrames, or (None, None) on failure.
    """
    add_rows, lai_rows = [], []

    for filename in sorted(os.listdir(hit_path)):
        if not filename.endswith('.glm.logistic.hybrid'):
            continue

        base = filename.split('.')[0]
        file_parts = base.split('_')
        window_start = int(file_parts[-2])
        window_end   = int(file_parts[-1])

        if allowed_windows is not None and (window_start, window_end) not in allowed_windows:
            continue

        filepath = os.path.join(hit_path, filename)
        try:
            df = pd.read_csv(filepath, sep='\t')
            df.columns = ['CHROM','POS','ID','REF','ALT','PROV_REF','A1','OMITTED',
                          'A1_FREQ','FIRTH','TEST','OBS_CT','OR','LOG_OR_SE',
                          'L95','U95','Z_STAT','P','ERRCODE']
        except Exception as e:
            print(f"  Error reading {filename}: {e}")
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


def run_analysis(all_add, all_lai, label):
    """
    Apply BY correction, merge ADD+LAI, compute stats.
    Returns (merged df with P_BY columns, stats dict).
    """
    _, add_by, _, _ = multipletests(all_add['P'].values, method='fdr_by')
    all_add = all_add.copy()
    all_add['P_BY'] = add_by

    _, lai_by, _, _ = multipletests(all_lai['P'].values, method='fdr_by')
    all_lai = all_lai.copy()
    all_lai['P_BY'] = lai_by

    lai_cols = all_lai[['POS','ID','window_start','P','P_BY']].rename(
        columns={'P': 'LAI_P', 'P_BY': 'LAI_P_BY'}
    )
    merged = pd.merge(all_add, lai_cols, on=['POS','ID','window_start'], how='inner')

    stats = {
        'n_snps':     len(merged),
        'n_add_sig':  (merged['P_BY'] < 0.05).sum(),
        'n_lai_sig':  (merged['LAI_P_BY'] < 0.05).sum(),
        'n_both':     ((merged['P_BY'] < 0.05) & (merged['LAI_P_BY'] < 0.05)).sum(),
        'n_none':     ((merged['P_BY'] >= 0.05) & (merged['LAI_P_BY'] >= 0.05)).sum(),
        'n_causal':   ((merged['P_BY'] < 0.05) & (merged['LAI_P_BY'] >= 0.05)).sum(),
        'label':      label,
    }
    return merged, stats


def print_stats(hit_dir, stats):
    print(f"\n{hit_dir} [{stats['label']}] — {stats['n_snps']} total SNPs:")
    print(f"  ADD significant (P_BY<0.05):           {stats['n_add_sig']}")
    print(f"  LAI significant (P_BY<0.05):           {stats['n_lai_sig']}")
    print(f"  Both significant:                      {stats['n_both']}")
    print(f"  Neither significant:                   {stats['n_none']}")
    print(f"  SNP sig + LAI not sig (candidates):    {stats['n_causal']}")


results = []

for hit_dir in sorted(os.listdir(fine_mapping_folder)):
    hit_path = os.path.join(fine_mapping_folder, hit_dir)
    if not os.path.isdir(hit_path):
        continue

    parts = hit_dir.split('_')
    hit_ancestry = parts[0]
    pheno        = parts[1]
    chr_num      = parts[2].replace('chr', '')

    # --- Step 1: analysis on significant windows only ---
    all_add, all_lai = load_hit_data(hit_path)
    if all_add is None:
        print(f"Skipping {hit_dir}: no data found")
        continue

    merged, stats = run_analysis(all_add, all_lai, label='sig windows only')
    print_stats(hit_dir, stats)

    candidates = merged[(merged['P_BY'] < 0.05) & (merged['LAI_P_BY'] >= 0.05)].copy()

    # --- Step 2: if no candidates, extend ±N_EXT windows ---
    if len(candidates) == 0:
        print(f"  -> No candidates in sig windows. Trying ±{N_EXT} extension windows...")

        hit_6wind_path = os.path.join(fine_mapping_6wind_folder, hit_dir)
        if not os.path.isdir(hit_6wind_path):
            print(f"  -> Extension folder not found: {hit_6wind_path}")
            continue

        # Get sorted list of extension windows (only those in 6wind, not in sig)
        ext_windows_sorted = []
        for fn in sorted(os.listdir(hit_6wind_path)):
            if not fn.endswith('.glm.logistic.hybrid'):
                continue
            fp = fn.split('.')[0].split('_')
            ext_windows_sorted.append((int(fp[-2]), int(fp[-1])))
        ext_windows_sorted.sort()

        # Sig window boundaries
        sig_windows = set()
        for fn in os.listdir(hit_path):
            if not fn.endswith('.glm.logistic.hybrid'):
                continue
            fp = fn.split('.')[0].split('_')
            sig_windows.add((int(fp[-2]), int(fp[-1])))

        min_sig = min(w[0] for w in sig_windows)
        max_sig = max(w[1] for w in sig_windows)

        # Take N_EXT upstream and N_EXT downstream from extension windows
        upstream   = [w for w in ext_windows_sorted if w[1] <= min_sig][-N_EXT:]
        downstream = [w for w in ext_windows_sorted if w[0] >= max_sig][:N_EXT]
        selected_ext = set(upstream) | set(downstream)

        print(f"  -> Using {len(upstream)} upstream + {len(sig_windows)} sig + {len(downstream)} downstream windows")

        # Load extension windows from 6wind folder, sig windows from sig folder
        all_add_ext, all_lai_ext = load_hit_data(hit_6wind_path, allowed_windows=selected_ext)
        if all_add_ext is None:
            all_add_ext = pd.DataFrame()
            all_lai_ext = pd.DataFrame()

        # Combine sig results (already loaded) with extension results
        all_add_ext = pd.concat([all_add, all_add_ext], ignore_index=True)
        all_lai_ext = pd.concat([all_lai, all_lai_ext], ignore_index=True)

        merged, stats = run_analysis(all_add_ext, all_lai_ext, label=f'sig ± {N_EXT} windows')
        print_stats(hit_dir, stats)

        candidates = merged[(merged['P_BY'] < 0.05) & (merged['LAI_P_BY'] >= 0.05)].copy()

    if len(candidates) == 0:
        print(f"  -> No causal candidates found")
        continue

    # Rank by |beta|
    candidates['beta'] = np.log(candidates['OR'])
    candidates = candidates.sort_values('beta', key=lambda x: x.abs(), ascending=False)

    top = candidates.iloc[0]
    print(f"  -> Top SNP: {top['ID']} chr{top['CHROM']}:{top['POS']} "
          f"OR={top['OR']:.3f} beta={np.log(top['OR']):.3f} "
          f"P_BY={top['P_BY']:.4g} LAI_P_BY={top['LAI_P_BY']:.4g}")

    results.append({
        'hit':           hit_dir,
        'ancestry':      hit_ancestry,
        'pheno':         pheno,
        'chr':           chr_num,
        'analysis':      stats['label'],
        'window_start':  top['window_start'],
        'window_end':    top['window_end'],
        'SNP_ID':        top['ID'],
        'POS':           top['POS'],
        'REF':           top['REF'],
        'ALT':           top['ALT'],
        'A1':            top['A1'],
        'OR':            round(top['OR'], 4),
        'beta':          round(top['beta'], 4),
        'SE':            round(top['LOG_OR_SE'], 4),
        'OR_L95':        round(top['L95'], 4),
        'OR_U95':        round(top['U95'], 4),
        'P_ADD':         top['P'],
        'P_ADD_BY':      round(top['P_BY'], 6),
        'LAI_P':         top['LAI_P'],
        'LAI_P_BY':      round(top['LAI_P_BY'], 4),
        'OBS_CT':        int(top['OBS_CT']),
    })

print(f"\n{'='*50}")
if results:
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved {len(results_df)} hits to {output_file}")
else:
    print("No causal candidates found.")
