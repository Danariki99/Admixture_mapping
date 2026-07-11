import os
import sys
import subprocess

# TEST version of borzoi_sed_windows.py — runs Borzoi SED (gene-level) for ALL
# SNPs in the fine-mapping windows (windows_all_hg38.vcf), across the 4 folds.
#
# NEEDS the SEPARATE Borzoi environment (TensorFlow) and, in practice, a GPU.
# Point BORZOI_DIR at the borzoi repo (with examples/saved_models + examples/hg38).
# Parameterized on result_folder.

if len(sys.argv) < 2:
    print("Usage: python borzoi_sed_windows_test.py <result_folder>")
    sys.exit(1)

result_folder = sys.argv[1]
BORZOI_DIR  = os.environ.get('BORZOI_DIR', '/private/home/rsmerigl/codes/cleaned_codes/borzoi')
OUTPUT_BASE = os.path.join(result_folder, 'borzoi_results', 'sed_windows')
VCF_FILE    = os.path.join(result_folder, 'borzoi_results', 'windows_all_hg38.vcf')

PARAMS_JSON  = os.path.join(BORZOI_DIR, 'examples/params_pred.json')
TARGETS_FILE = os.path.join(BORZOI_DIR, 'examples/targets_human.txt')
BORZOI_SED   = os.path.join(BORZOI_DIR, 'src/scripts/borzoi_sed.py')
FOLDS = ['f3c0', 'f3c1', 'f3c2', 'f3c3']

if not os.path.exists(VCF_FILE):
    raise FileNotFoundError(f'VCF not found: {VCF_FILE}. Run windows_snps_vcf_creation_test.py first.')

os.makedirs(OUTPUT_BASE, exist_ok=True)
env = os.environ.copy()
env['BORZOI_HG38'] = os.path.join(BORZOI_DIR, 'examples/hg38')
env['BORZOI_DIR']  = BORZOI_DIR

for fold in FOLDS:
    model_h5   = os.path.join(BORZOI_DIR, f'examples/saved_models/{fold}/train/model0_best.h5')
    output_dir = os.path.join(OUTPUT_BASE, fold)
    os.makedirs(output_dir, exist_ok=True)
    if os.path.exists(os.path.join(output_dir, 'sed.h5')):
        print(f'[{fold}] sed.h5 already exists, skipping')
        continue
    cmd = ['python', BORZOI_SED, '--rc', '--stats', 'SED', '-t', TARGETS_FILE,
           '-o', output_dir, PARAMS_JSON, model_h5, VCF_FILE]
    print(f'\n[{fold}] Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, check=False, env=env)
    print(f'[{fold}] {"Done" if result.returncode == 0 else f"ERROR {result.returncode}"} -> {output_dir}')

print(f'\nAll folds complete. Output: {OUTPUT_BASE}/')
