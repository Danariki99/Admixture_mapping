"""
Runs Borzoi SED (gene-level) for ALL SNPs in the significant fine-mapping windows
(windows_all_hg38.vcf, ~8k SNPs), across the 4 model replicates.

Extends the candidate-only runs: covers every SNP tested in the significant
windows, not just the fine-mapping candidates.

Output: borzoi_results/sed_windows/{f3c0..f3c3}/sed.h5
"""

import os
import subprocess

BORZOI_DIR  = '/private/home/rsmerigl/codes/cleaned_codes/borzoi'
OUTPUT_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/sed_windows'
VCF_FILE    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/windows_all_hg38.vcf'

PARAMS_JSON  = os.path.join(BORZOI_DIR, 'examples/params_pred.json')
TARGETS_FILE = os.path.join(BORZOI_DIR, 'examples/targets_human.txt')
BORZOI_SED   = os.path.join(BORZOI_DIR, 'src/scripts/borzoi_sed.py')

FOLDS = ['f3c0', 'f3c1', 'f3c2', 'f3c3']

# When run as a SLURM array (sbatch --array=0-3), each task does ONE fold, so
# every fold gets the full time budget and none is truncated. Otherwise all 4
# folds run sequentially in this single process.
_aid = os.environ.get('SLURM_ARRAY_TASK_ID')
if _aid is not None:
    FOLDS = [FOLDS[int(_aid)]]

os.makedirs(OUTPUT_BASE, exist_ok=True)

if not os.path.exists(VCF_FILE):
    raise FileNotFoundError(f'VCF not found: {VCF_FILE}. Run windows_snps_vcf_creation.py first.')

print(f'Using VCF: {VCF_FILE}')

env = os.environ.copy()
env['BORZOI_HG38'] = os.path.join(BORZOI_DIR, 'examples/hg38')
env['BORZOI_DIR']  = BORZOI_DIR

for fold in FOLDS:
    model_h5   = os.path.join(BORZOI_DIR, f'examples/saved_models/{fold}/train/model0_best.h5')
    output_dir = os.path.join(OUTPUT_BASE, fold)
    os.makedirs(output_dir, exist_ok=True)

    sed_out = os.path.join(output_dir, 'sed.h5')
    if os.path.exists(sed_out):
        print(f'[{fold}] sed.h5 already exists, skipping')
        continue

    cmd = [
        'python', BORZOI_SED,
        '--rc',
        '--stats', 'SED',
        '-t', TARGETS_FILE,
        '-o', output_dir,
        PARAMS_JSON,
        model_h5,
        VCF_FILE,
    ]

    print(f'\n[{fold}] Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, check=False, env=env)
    if result.returncode != 0:
        print(f'[{fold}] ERROR: exit code {result.returncode}')
    else:
        print(f'[{fold}] Done → {output_dir}')

print('\nAll folds complete.')
print(f'Output: {OUTPUT_BASE}/')
