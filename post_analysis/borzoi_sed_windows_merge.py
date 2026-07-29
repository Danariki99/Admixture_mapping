"""
Merge the sharded Borzoi SED outputs back into one sed.h5 per fold.

borzoi_sed_windows.py (SLURM array) writes sed_windows/{fold}/shard{NN}/sed.h5,
each covering a contiguous slice of the VCF. This concatenates the shards of a
fold into sed_windows/{fold}/sed.h5 with the SAME layout as an unsharded run, so
borzoi_sed_post_processing.py reads it unchanged.

Layout of each sed.h5 (two axes):
  per-SNP  (n_snp):   snp, chr, pos, ref_allele, alt_allele
  per-pair (n_pairs): SED (n_pairs, n_tracks), gene, si  (si -> snp index)
  global:             target_ids, target_labels
Merging concatenates both axes and offsets each shard's `si` by the cumulative
SNP count so the pair->SNP mapping stays correct.
"""

import os
import glob
import h5py
import numpy as np

SED_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/sed_windows'
FOLDS    = ['f3c0', 'f3c1', 'f3c2', 'f3c3']

PER_SNP  = ['snp', 'chr', 'pos', 'ref_allele', 'alt_allele']
PER_PAIR = ['SED', 'gene']          # 'si' handled separately (needs offset)
GLOBAL   = ['target_ids', 'target_labels']


def merge_fold(fold):
    fold_dir = os.path.join(SED_BASE, fold)
    shards = sorted(glob.glob(os.path.join(fold_dir, 'shard*', 'sed.h5')))
    if not shards:
        print(f'[{fold}] no shards found, skipping')
        return

    per_snp  = {k: [] for k in PER_SNP}
    per_pair = {k: [] for k in PER_PAIR}
    si_parts = []
    global_data = {}
    snp_offset = 0

    for sp in shards:
        with h5py.File(sp, 'r') as f:
            for k in PER_SNP:
                per_snp[k].append(f[k][()])
            for k in PER_PAIR:
                per_pair[k].append(f[k][()])
            si_parts.append(f['si'][()] + snp_offset)   # offset into the merged SNP axis
            snp_offset += f['snp'].shape[0]
            if not global_data:
                for k in GLOBAL:
                    global_data[k] = f[k][()]

    out_path = os.path.join(fold_dir, 'sed.h5')
    with h5py.File(out_path, 'w') as out:
        for k in PER_SNP:
            out.create_dataset(k, data=np.concatenate(per_snp[k]))
        for k in PER_PAIR:
            out.create_dataset(k, data=np.concatenate(per_pair[k], axis=0))
        out.create_dataset('si', data=np.concatenate(si_parts))
        for k in GLOBAL:
            out.create_dataset(k, data=global_data[k])

    with h5py.File(out_path, 'r') as chk:
        print(f'[{fold}] merged {len(shards)} shards -> {out_path}  '
              f'(n_snp={chk["snp"].shape[0]}, n_pairs={chk["SED"].shape[0]}, '
              f'max_si={int(chk["si"][()].max())})')


if __name__ == '__main__':
    for fold in FOLDS:
        merge_fold(fold)
    print('\nDone. sed_windows/{fold}/sed.h5 ready for borzoi_sed_post_processing.py')
