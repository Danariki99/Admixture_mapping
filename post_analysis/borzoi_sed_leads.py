"""
Runs Borzoi SED (SNP Expression Difference, gene-level) for the 5 lead/causal
candidate SNPs, across the 4 model replicates.

Unlike SAD (track-level), SED uses the gene GTF to compute a per-GENE predicted
expression change for every gene whose exons fall in the prediction window of
each SNP — i.e. "which genes in the region, and how their expression changes".

Output: borzoi_results/sed_leads/{f3c0..f3c3}/sed.h5
"""

import os
import subprocess

BORZOI_DIR  = '/private/home/rsmerigl/codes/cleaned_codes/borzoi'
OUTPUT_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/sed_leads'
VCF_FILE    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/leads_5_hg38.vcf'

PARAMS_JSON  = os.path.join(BORZOI_DIR, 'examples/params_pred.json')
TARGETS_FILE = os.path.join(BORZOI_DIR, 'examples/targets_human.txt')
BORZOI_SED   = os.path.join(BORZOI_DIR, 'src/scripts/borzoi_sed.py')

FOLDS = ['f3c0', 'f3c1', 'f3c2', 'f3c3']

os.makedirs(OUTPUT_BASE, exist_ok=True)

if not os.path.exists(VCF_FILE):
    raise FileNotFoundError(f'VCF not found: {VCF_FILE}')

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
