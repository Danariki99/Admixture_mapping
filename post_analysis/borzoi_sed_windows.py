"""
Runs Borzoi SED (gene-level) for ALL SNPs in the significant fine-mapping windows
(windows_all_hg38.vcf, ~8k SNPs), across the 4 model replicates.

Extends the candidate-only runs: covers every SNP tested in the significant
windows, not just the fine-mapping candidates.

SHARDED execution (SLURM array): running all ~8k SNPs of a fold in one task is
too slow for a single wall-time budget, so the work is split over
FOLD x SHARD tasks. With --array=0-(4*N_SHARDS-1):
    fold_idx  = SLURM_ARRAY_TASK_ID // N_SHARDS
    shard_idx = SLURM_ARRAY_TASK_ID %  N_SHARDS
Each task scores one contiguous slice of the (already sorted, hg38-matched) VCF
and writes sed_windows/{fold}/shard{NN}/sed.h5. Merge them afterwards with
borzoi_sed_windows_merge.py -> sed_windows/{fold}/sed.h5.

N_SHARDS defaults to 20 (~400 SNPs/task) and can be overridden with the
BORZOI_N_SHARDS env var; it MUST match the --array range and the merge script.

Without SLURM_ARRAY_TASK_ID this falls back to the old behaviour: all 4 folds,
all SNPs, sequentially in one process.

Output: borzoi_results/sed_windows/{f3c0..f3c3}/shard{NN}/sed.h5
"""

import os
import subprocess

BORZOI_DIR  = '/private/home/rsmerigl/codes/cleaned_codes/borzoi'
OUTPUT_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/sed_windows'
VCF_FILE    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/windows_all_hg38.vcf'

PARAMS_JSON  = os.path.join(BORZOI_DIR, 'examples/params_pred.json')
TARGETS_FILE = os.path.join(BORZOI_DIR, 'examples/targets_human.txt')
BORZOI_SED   = os.path.join(BORZOI_DIR, 'src/scripts/borzoi_sed.py')

ALL_FOLDS = ['f3c0', 'f3c1', 'f3c2', 'f3c3']
N_SHARDS  = int(os.environ.get('BORZOI_N_SHARDS', '20'))

if not os.path.exists(VCF_FILE):
    raise FileNotFoundError(f'VCF not found: {VCF_FILE}. Run windows_snps_vcf_creation.py first.')

os.makedirs(OUTPUT_BASE, exist_ok=True)


def read_vcf(path):
    """Return (header_lines, snp_lines) preserving order."""
    header, snps = [], []
    with open(path) as fh:
        for line in fh:
            (header if line.startswith('#') else snps).append(line)
    return header, snps


def contiguous_slice(n, n_shards, k):
    """The k-th of n_shards contiguous, as-even-as-possible index ranges of n items."""
    base, rem = divmod(n, n_shards)
    start = k * base + min(k, rem)
    stop  = start + base + (1 if k < rem else 0)
    return start, stop


def run_borzoi(fold, output_dir, vcf_path):
    model_h5 = os.path.join(BORZOI_DIR, f'examples/saved_models/{fold}/train/model0_best.h5')
    os.makedirs(output_dir, exist_ok=True)
    sed_out = os.path.join(output_dir, 'sed.h5')
    if os.path.exists(sed_out):
        print(f'[{fold}] {output_dir}: sed.h5 already exists, skipping')
        return
    cmd = ['python', BORZOI_SED, '--rc', '--stats', 'SED',
           '-t', TARGETS_FILE, '-o', output_dir,
           PARAMS_JSON, model_h5, vcf_path]
    print(f'\n[{fold}] Running: {" ".join(cmd)}')
    env = os.environ.copy()
    env['BORZOI_HG38'] = os.path.join(BORZOI_DIR, 'examples/hg38')
    env['BORZOI_DIR']  = BORZOI_DIR
    result = subprocess.run(cmd, check=False, env=env)
    if result.returncode != 0:
        print(f'[{fold}] ERROR: exit code {result.returncode}')
    else:
        print(f'[{fold}] Done -> {output_dir}')


header, snps = read_vcf(VCF_FILE)
print(f'Using VCF: {VCF_FILE}  ({len(snps)} SNPs)')

_aid = os.environ.get('SLURM_ARRAY_TASK_ID')

if _aid is not None:
    # ── sharded array mode ────────────────────────────────────────────────────
    aid = int(_aid)
    fold_idx, shard_idx = divmod(aid, N_SHARDS)
    if fold_idx >= len(ALL_FOLDS):
        raise ValueError(f'array id {aid} -> fold {fold_idx} out of range '
                         f'(N_SHARDS={N_SHARDS}, expected --array=0-{len(ALL_FOLDS)*N_SHARDS-1})')
    fold = ALL_FOLDS[fold_idx]
    lo, hi = contiguous_slice(len(snps), N_SHARDS, shard_idx)
    shard_snps = snps[lo:hi]
    output_dir = os.path.join(OUTPUT_BASE, fold, f'shard{shard_idx:02d}')
    os.makedirs(output_dir, exist_ok=True)
    shard_vcf = os.path.join(output_dir, 'variants.vcf')
    with open(shard_vcf, 'w') as out:
        out.writelines(header)
        out.writelines(shard_snps)
    print(f'array {aid} -> fold {fold}  shard {shard_idx}/{N_SHARDS}  '
          f'SNPs [{lo}:{hi}] (n={len(shard_snps)})')
    run_borzoi(fold, output_dir, shard_vcf)
    print(f'\nShard done. Merge with borzoi_sed_windows_merge.py when all tasks finish.')
else:
    # ── fallback: all folds, all SNPs, sequential ─────────────────────────────
    for fold in ALL_FOLDS:
        run_borzoi(fold, os.path.join(OUTPUT_BASE, fold), VCF_FILE)
    print('\nAll folds complete.')
    print(f'Output: {OUTPUT_BASE}/')
