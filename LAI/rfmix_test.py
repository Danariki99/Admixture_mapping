"""
TEST: Local Ancestry Inference with RFMix v2 (replaces the GNomix step).

For each per-chromosome query VCF (produced by chrom_division.py in
result_folder/vcf_folder/chrN.vcf.gz) it runs the RFMix binary bundled OUTSIDE
Admixture_mapping at cleaned_codes/rfmix/rfmix (built from slowkoni/rfmix, like
gnomix/ and plink2), using:
  reference   = data_folder/panel_chrom/chrN.vcf.gz
  sample map  = data_folder/sample_map_filtered.txt   (sample <TAB> ancestry)
  genetic map = all data_folder/maps/chrN.gmap concatenated (chr pos cM)

RFMix writes chrN.msp.tsv (+ .fb.tsv, .Q, .sis.tsv); we rename chrN.msp.tsv ->
chrN.msp so the rest of the pipeline finds it exactly like the GNomix output
(this replaces LAI/files_moving.sh).

Same interface as gnomix_training_test.py:
  --vcf_folder <result_folder/vcf_folder>  --data_folder  --result_folder
"""

import os
import re
import sys
import glob
import argparse
import subprocess

# RFMix lives in cleaned_codes/rfmix/ (sibling of Admixture_mapping), like gnomix/
RFMIX = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'rfmix', 'rfmix'))


def chr_key(path):
    m = re.search(r'chr(\w+)', os.path.basename(path))
    tok = m.group(1) if m else ''
    return (0, int(tok)) if tok.isdigit() else (1, tok)


def build_genetic_map(maps_dir, out_gmap):
    """Concatenate the per-chromosome gmaps (already chr/pos/cM) into one
       whole-genome genetic map, as RFMix expects."""
    files = sorted(glob.glob(os.path.join(maps_dir, 'chr*.gmap')), key=chr_key)
    with open(out_gmap, 'w') as o:
        for f in files:
            with open(f) as fi:
                for line in fi:
                    if line.strip() and not line.startswith('#'):
                        o.write(line)
    return out_gmap


def ensure_index(vcf):
    if not (os.path.exists(vcf + '.tbi') or os.path.exists(vcf + '.csi')):
        subprocess.run(['tabix', '-p', 'vcf', vcf], check=True)


def main():
    parser = argparse.ArgumentParser(description='Run RFMix per chromosome on the test data')
    parser.add_argument('--vcf_folder', required=True, help='per-chromosome query VCFs (chrN.vcf.gz)')
    parser.add_argument('--data_folder', required=True)
    parser.add_argument('--result_folder', required=True)
    args = parser.parse_args()

    if not os.path.exists(RFMIX):
        sys.exit(f'RFMix binary not found: {RFMIX}\n  Build it: see the README (clone slowkoni/rfmix into cleaned_codes/rfmix and make).')

    maps_dir   = os.path.join(args.data_folder, 'maps')
    panel_dir  = os.path.join(args.data_folder, 'panel_chrom')
    sample_map = os.path.join(args.data_folder, 'sample_map_filtered.txt')
    msp_folder = os.path.join(args.result_folder, 'msp_folder')
    tmp_folder = os.path.join(args.result_folder, 'tmp')
    os.makedirs(msp_folder, exist_ok=True)
    os.makedirs(tmp_folder, exist_ok=True)

    gmap = build_genetic_map(maps_dir, os.path.join(tmp_folder, 'rfmix_allchrs.gmap'))
    print(f'Genetic map (whole genome): {gmap}')

    queries = sorted(glob.glob(os.path.join(args.vcf_folder, 'chr*.vcf.gz')), key=chr_key)
    for query in queries:
        chrom = os.path.basename(query).replace('.vcf.gz', '')          # e.g. chr6
        ref = os.path.join(panel_dir, f'{chrom}.vcf.gz')
        if not os.path.exists(ref):
            print(f'  {chrom}: reference panel {ref} missing, skipping'); continue

        ensure_index(query)
        ensure_index(ref)

        out_prefix = os.path.join(msp_folder, chrom)
        cmd = [RFMIX, '-f', query, '-r', ref, '-m', sample_map, '-g', gmap,
               '-o', out_prefix, '--n-threads=8', f'--chromosome={chrom}']
        print(f'\n[{chrom}] {" ".join(cmd)}')
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            print(f'  {chrom}: RFMix FAILED'); continue

        # uniform output name for the rest of the pipeline: chrN.msp.tsv -> chrN.msp
        msp_tsv = out_prefix + '.msp.tsv'
        if os.path.exists(msp_tsv):
            os.replace(msp_tsv, out_prefix + '.msp')
            print(f'  -> {out_prefix}.msp')
        else:
            print(f'  {chrom}: expected {msp_tsv} not produced')


if __name__ == '__main__':
    main()
