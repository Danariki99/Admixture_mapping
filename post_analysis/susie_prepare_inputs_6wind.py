"""
Prepares SuSiE inputs for the ±6-window extended analysis (UKBB in-sample LD).
Fix 2: monomorphic/NaN SNPs are removed (not patched).
Fix 3: strand-ambiguous SNPs (A/T, C/G) are removed; REF/ALT saved in output.
Z-scores are loaded from both fine_mapping_new/ and fine_mapping_new_6wind/,
then combined and deduplicated by best p-value.
"""

import os
import numpy as np
import pandas as pd

fine_mapping_folder    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
fine_mapping_6w_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new_6wind'
susie_ld_folder        = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld_6wind'
output_folder          = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_inputs_6wind'

os.makedirs(output_folder, exist_ok=True)


def is_strand_ambiguous(ref, alt):
    return frozenset({ref.upper(), alt.upper()}) in (frozenset({'A', 'T'}), frozenset({'C', 'G'}))


def load_zscores(hit_label):
    rows = []
    for folder in [fine_mapping_folder, fine_mapping_6w_folder]:
        hit_path = os.path.join(folder, hit_label)
        if not os.path.isdir(hit_path):
            continue
        for filename in sorted(os.listdir(hit_path)):
            if not filename.endswith('.glm.logistic.hybrid'):
                continue
            filepath = os.path.join(hit_path, filename)
            try:
                df = pd.read_csv(filepath, sep='\t')
                df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'PROV_REF', 'A1', 'OMITTED',
                              'A1_FREQ', 'FIRTH', 'TEST', 'OBS_CT', 'OR', 'LOG_OR_SE',
                              'L95', 'U95', 'Z_STAT', 'P', 'ERRCODE']
            except Exception as e:
                print(f"  Error reading {filename}: {e}")
                continue
            add = df[df['TEST'] == 'ADD'].copy()
            for col in ['Z_STAT', 'P']:
                add[col] = pd.to_numeric(add[col], errors='coerce')
            add = add.dropna(subset=['Z_STAT', 'P'])
            rows.append(add[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Z_STAT', 'P', 'OBS_CT']])

    if not rows:
        return None

    all_snps = pd.concat(rows, ignore_index=True)
    all_snps = (all_snps
                .sort_values('P')
                .drop_duplicates(subset='ID', keep='first')
                .sort_values(['CHROM', 'POS'])
                .reset_index(drop=True))
    all_snps = all_snps.rename(columns={'Z_STAT': 'Z', 'OBS_CT': 'N'})
    return all_snps[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Z', 'N']]


def load_ld_plink(vcor1_file, vars_file):
    snp_ids = pd.read_csv(vars_file, header=None)[0].tolist()
    R = np.loadtxt(vcor1_file)
    return snp_ids, R


for hit_dir in sorted(os.listdir(fine_mapping_folder)):
    hit_path = os.path.join(fine_mapping_folder, hit_dir)
    if not os.path.isdir(hit_path):
        continue

    hit_label  = hit_dir
    vcor1_file = os.path.join(susie_ld_folder, hit_label, f'{hit_label}_ld_6wind.phased.vcor1')
    vars_file  = os.path.join(susie_ld_folder, hit_label, f'{hit_label}_ld_6wind.phased.vcor1.vars')

    if not os.path.exists(vcor1_file):
        print(f"{hit_label}: 6wind LD file not found, skipping")
        continue

    print(f"\n{hit_label}")

    zscores = load_zscores(hit_label)
    if zscores is None:
        print(f"  No z-scores found, skipping")
        continue
    print(f"  Z-scores: {len(zscores)} SNPs (orig + 6wind combined)")

    ld_snp_ids, R = load_ld_plink(vcor1_file, vars_file)
    print(f"  LD matrix: {R.shape[0]} SNPs")

    # Align by ID (same UKBB VCF for both GWAS and LD)
    shared      = set(zscores['ID']) & set(ld_snp_ids)
    if len(shared) == 0:
        print(f"  No shared SNPs, skipping")
        continue

    zscores_aln = zscores[zscores['ID'].isin(shared)].reset_index(drop=True)
    id_to_idx   = {sid: i for i, sid in enumerate(ld_snp_ids)}
    col_idx     = [id_to_idx[sid] for sid in zscores_aln['ID']]
    R_aln       = R[np.ix_(col_idx, col_idx)]
    print(f"  Shared SNPs: {len(zscores_aln)}")

    # Fix 3: remove strand-ambiguous SNPs (A/T or C/G — orientation unresolvable)
    sa_mask = zscores_aln.apply(lambda r: is_strand_ambiguous(r['REF'], r['ALT']), axis=1).values
    n_ambig = int(sa_mask.sum())
    if n_ambig > 0:
        print(f"  Strand-ambiguous SNPs removed: {n_ambig}")
        keep        = ~sa_mask
        zscores_aln = zscores_aln[keep].reset_index(drop=True)
        R_aln       = R_aln[np.ix_(keep, keep)]

    # Fix 2: symmetrize, then remove monomorphic and partial-NaN SNPs
    R_aln        = (R_aln + R_aln.T) / 2
    diag         = np.diag(R_aln)
    mono_mask    = np.isnan(diag)
    partial_mask = (~mono_mask) & np.isnan(R_aln).any(axis=1)
    n_mono       = int(mono_mask.sum())
    n_partial    = int(partial_mask.sum())
    if n_mono > 0:
        print(f"  Removed monomorphic SNPs (all-NaN in LD): {n_mono}")
    if n_partial > 0:
        print(f"  Removed SNPs with partial NaN in LD: {n_partial}")
    remove_mask = mono_mask | partial_mask
    if remove_mask.any():
        keep        = ~remove_mask
        zscores_aln = zscores_aln[keep].reset_index(drop=True)
        R_aln       = R_aln[np.ix_(keep, keep)]

    if len(zscores_aln) == 0:
        print(f"  No SNPs remaining after filtering, skipping")
        continue

    np.fill_diagonal(R_aln, 1.0)
    print(f"  Final SNPs: {len(zscores_aln)}")

    hit_out = os.path.join(output_folder, hit_label)
    os.makedirs(hit_out, exist_ok=True)

    zscores_aln.to_csv(os.path.join(hit_out, f'{hit_label}_zscores.tsv'), sep='\t', index=False)
    np.savetxt(os.path.join(hit_out, f'{hit_label}_ld.gz'), R_aln, fmt='%.6f')
    with open(os.path.join(hit_out, f'{hit_label}_snp_ids.txt'), 'w') as f:
        f.write('\n'.join(zscores_aln['ID'].tolist()) + '\n')

    print(f"  Saved → {hit_out}/")

print("\nDone.")
