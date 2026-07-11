"""
Entropy "mountain" plot per hit (EXTENDED, ±6 windows):
  - top   : global (all-samples) LD r^2 triangle  [SuSiE_ld_6wind]
  - bottom: per-locus HAPLOTYPE entropy, computed on all KING-cutoff samples.

Entropy of each locus = Shannon entropy (bits) of the distinct phased haplotypes
observed in the ±7-SNP window centred on that locus (15 SNPs), over all
2 × 415,792 haplotypes.

Haplotypes are extracted on the fly from the phased UKBB VCF:
  bcftools view -r <region>          (fast, tabix-indexed slice)
  plink2 --keep <kcutoff> --extract <hit .vars> --export haps

Both the LD triangle and the entropy are on the SAME all-samples cohort
(LD globale + entropia globale).

Output: ld_mountain_plots_all/{hit}__allsamples__entropy.png / .pdf
"""

import os
import glob
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

BASE     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb'
LD_DIR   = f'{BASE}/SuSiE_ld_6wind'
FM_DIRS  = [f'{BASE}/fine_mapping_new', f'{BASE}/fine_mapping_new_6wind']
VCF      = f'{BASE}/vcf_file/ukbb.vcf.gz'
KCUTOFF  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
PLINK2   = '/private/home/rsmerigl/plink2'
OUT_DIR  = os.path.join(os.path.dirname(__file__), 'ld_mountain_plots_all')
WORK     = os.path.join(OUT_DIR, '_entropy_work')
# shared, persistent entropy cache (same one used by ld_combined_plots.py)
ENTROPY_CACHE = f'{BASE}/entropy_cache'

HALF = 7        # ±7 SNP window (15 SNPs)

mpl.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans']})


def fm_positions(hit):
    """POS per rsID from fine-mapping ADD rows (for the region span + x labels)."""
    pos = {}
    for fm_dir in FM_DIRS:
        for f in glob.glob(os.path.join(fm_dir, hit, '*.glm.logistic.hybrid')):
            d = pd.read_csv(f, sep='\t', usecols=['ID', 'POS', 'TEST'])
            d = d[d['TEST'] == 'ADD']
            pos.update(dict(zip(d['ID'], pd.to_numeric(d['POS'], errors='coerce'))))
    return pos


def extract_haplotypes(hit, chrom, region, var_file, tag):
    """bcftools slice -> plink2 haps export. Returns path to the .haps file."""
    os.makedirs(WORK, exist_ok=True)
    region_vcf = os.path.join(WORK, f'{tag}_region.vcf.gz')
    out_base   = os.path.join(WORK, tag)
    subprocess.run(['bcftools', 'view', '-r', region, VCF, '-Oz', '-o', region_vcf], check=True)
    subprocess.run(['tabix', '-f', '-p', 'vcf', region_vcf], check=True)
    subprocess.run([PLINK2, '--vcf', region_vcf, '--keep', KCUTOFF,
                    '--extract', var_file, '--export', 'haps', '--out', out_base],
                   check=True, stdout=subprocess.DEVNULL)
    return out_base + '.haps', region_vcf, out_base


def parse_haps(haps_file):
    """Return (ids, H) where H is a variants × (2N) uint8 allele matrix."""
    ids, rows = [], []
    with open(haps_file) as f:
        for line in f:
            p = line.split(' ', 5)
            ids.append(p[1])
            a = np.frombuffer(p[5].rstrip('\n').encode(), np.uint8)[::2] - 48
            rows.append(a.astype(np.uint8))
    return ids, np.vstack(rows)


def haplotype_entropy(Hal, present):
    """Per-locus Shannon entropy (bits) over the ±HALF window of present SNPs."""
    V, twoN = Hal.shape
    ent = np.full(V, np.nan)
    for i in range(V):
        idx = [j for j in range(max(0, i - HALF), min(V, i + HALF + 1)) if present[j]]
        if len(idx) < 2:
            continue
        sub = Hal[idx]
        codes = (1 << np.arange(len(idx))).astype(np.int64) @ sub
        _, c = np.unique(codes, return_counts=True)
        p = c / twoN
        ent[i] = -(p * np.log2(p)).sum()
    return ent


def make_plot(hit):
    ld_file  = os.path.join(LD_DIR, hit, f'{hit}_ld_6wind.phased.vcor1')
    var_file = ld_file + '.vars'
    if not os.path.exists(ld_file):
        print(f'  {hit}: no LD matrix, skip'); return

    snps = pd.read_csv(var_file, header=None)[0].astype(str).tolist()
    N = len(snps)
    posmap = fm_positions(hit)
    pos = np.array([posmap.get(s, np.nan) for s in snps], dtype=float)
    finite_pos = pos[np.isfinite(pos)]
    if finite_pos.size == 0:
        print(f'  {hit}: no positions, skip'); return
    chrom = hit.split('_')[-1].replace('chr', '')
    region = f'{chrom}:{int(finite_pos.min())}-{int(finite_pos.max())}'

    cf = os.path.join(ENTROPY_CACHE, f'{hit}_entropy.tsv')
    if os.path.exists(cf):
        print(f'  {hit}: {N} SNPs — reusing entropy cache {os.path.basename(cf)}')
        ent = pd.read_csv(cf, sep='\t').set_index('ID')['entropy'].reindex(snps).to_numpy(dtype=float)
    else:
        print(f'  {hit}: {N} SNPs — extracting haplotypes ({region}) ...')
        haps_file = region_vcf = out_base = None
        try:
            haps_file, region_vcf, out_base = extract_haplotypes(hit, chrom, region, var_file, hit)
            ids, H = parse_haps(haps_file)
            print(f'    haps {H.shape}, computing entropy ...')
            idpos = {v: i for i, v in enumerate(ids)}
            order = [idpos.get(v, -1) for v in snps]
            present = np.array([o >= 0 for o in order])
            Hal = np.vstack([H[o] if o >= 0 else np.zeros(H.shape[1], np.uint8) for o in order])
            ent = haplotype_entropy(Hal, present)
        finally:
            # clean up the big intermediates for this hit
            for pth in (region_vcf, (region_vcf + '.tbi') if region_vcf else None,
                        (out_base + '.haps') if out_base else None,
                        (out_base + '.sample') if out_base else None,
                        (out_base + '.log') if out_base else None):
                if pth and os.path.exists(pth):
                    os.remove(pth)
        os.makedirs(ENTROPY_CACHE, exist_ok=True)
        pd.DataFrame({'ID': snps, 'entropy': ent}).to_csv(cf, sep='\t', index=False)
        print(f'    entropy cached -> {os.path.basename(cf)}')

    R = np.loadtxt(ld_file)
    ld2 = R ** 2

    # ── figure ──────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(9, 5.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.4], width_ratios=[45, 1],
                          hspace=0.05, wspace=0.02)
    axm = fig.add_subplot(gs[0, 0])
    axp = fig.add_subplot(gs[1, 0], sharex=axm)
    cax = fig.add_subplot(gs[0, 1])

    ii, jj = np.meshgrid(np.arange(N + 1), np.arange(N + 1), indexing='ij')
    X = (ii + jj) / 2.0
    Y = (jj - ii) / 2.0
    C = ld2.copy().astype(float)
    C[np.tril_indices(N, -1)] = np.nan
    pc = axm.pcolormesh(X, Y, C, cmap='Reds', vmin=0, vmax=1, shading='flat', rasterized=True)
    axm.set_ylim(0, N / 2.0); axm.set_xlim(0, N); axm.set_yticks([])
    axm.tick_params(axis='x', labelbottom=False, bottom=False)
    for s in ('top', 'right', 'left'):
        axm.spines[s].set_visible(False)
    axm.set_title(f"{hit.replace('NAT', 'AMR')}  (all samples LD + haplotype entropy, ±6 windows)",
                  fontsize=9, fontweight='bold')
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r'LD $r^2$', fontsize=7); cb.ax.tick_params(labelsize=6)

    x = np.arange(N)
    axp.fill_between(x, 0, np.nan_to_num(ent, nan=0.0), where=np.isfinite(ent),
                     color='#6a3d9a', alpha=0.25, linewidth=0)
    axp.plot(x, ent, color='#6a3d9a', lw=0.6, rasterized=True)
    axp.set_ylabel('Haplotype entropy\n(±7 SNP, bits)', fontsize=8)
    axp.set_ylim(bottom=0)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)

    nt = 6
    idx = np.linspace(0, N - 1, nt).astype(int)
    labs = [f'{pos[i]/1e6:.2f}' if np.isfinite(pos[i]) else '' for i in idx]
    axp.set_xticks(idx); axp.set_xticklabels(labs, fontsize=6)
    axp.set_xlabel(f'Position on chr{chrom} (Mb)', fontsize=8)
    axp.tick_params(axis='y', labelsize=6)

    os.makedirs(OUT_DIR, exist_ok=True)
    base = os.path.join(OUT_DIR, f'{hit}__allsamples__entropy')
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{base}.pdf', dpi=300, bbox_inches='tight', format='pdf')
    plt.close(fig)
    print(f'    saved -> {os.path.basename(base)}.png/.pdf')


if __name__ == '__main__':
    import sys
    hits = sys.argv[1:] or sorted(
        d for d in os.listdir(LD_DIR) if os.path.isdir(os.path.join(LD_DIR, d)))
    for hit in hits:
        try:
            make_plot(hit)
        except Exception as e:
            print(f'  {hit}: FAILED ({e})')
    if os.path.isdir(WORK) and not os.listdir(WORK):
        shutil.rmtree(WORK, ignore_errors=True)
