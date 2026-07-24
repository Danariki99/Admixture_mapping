"""
TEST run of Borzoi SED on a SMALL, SORTED VCF (borzoi_results/sed_test/), ONE fold
only, to validate that the SNPs which came out all-zero in the full sed_windows
run now produce non-zero SED — before committing to the multi-day full run.

Identical Borzoi call as borzoi_sed_windows.py; only the VCF (5 SNPs, sorted) and
the fold list differ.
"""

import os
import subprocess

BORZOI_DIR  = '/private/home/rsmerigl/codes/cleaned_codes/borzoi'
OUTPUT_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/sed_test'
VCF_FILE    = os.path.join(OUTPUT_BASE, 'test_snps_hg38.vcf')

PARAMS_JSON  = os.path.join(BORZOI_DIR, 'examples/params_pred.json')
TARGETS_FILE = os.path.join(BORZOI_DIR, 'examples/targets_human.txt')
BORZOI_SED   = os.path.join(BORZOI_DIR, 'src/scripts/borzoi_sed.py')

FOLDS = ['f3c0']         # just one fold for the test

os.makedirs(OUTPUT_BASE, exist_ok=True)
if not os.path.exists(VCF_FILE):
    raise FileNotFoundError(f'Test VCF not found: {VCF_FILE}')
print(f'Using TEST VCF: {VCF_FILE}')

env = os.environ.copy()
env['BORZOI_HG38'] = os.path.join(BORZOI_DIR, 'examples/hg38')
env['BORZOI_DIR']  = BORZOI_DIR

for fold in FOLDS:
    model_h5   = os.path.join(BORZOI_DIR, f'examples/saved_models/{fold}/train/model0_best.h5')
    output_dir = os.path.join(OUTPUT_BASE, fold)
    os.makedirs(output_dir, exist_ok=True)
    sed_out = os.path.join(output_dir, 'sed.h5')
    if os.path.exists(sed_out):
        print(f'[{fold}] sed.h5 exists, removing for a clean test'); os.remove(sed_out)

    cmd = ['python', BORZOI_SED, '--rc', '--stats', 'SED', '-t', TARGETS_FILE,
           '-o', output_dir, PARAMS_JSON, model_h5, VCF_FILE]
    print(f'\n[{fold}] Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, check=False, env=env)
    print(f'[{fold}] {"Done" if result.returncode == 0 else f"ERROR {result.returncode}"} -> {output_dir}')

print('\nTest complete. Verify sed_test/f3c0/sed.h5 has non-zero SED for rs3177928, rs6906021, etc.')
