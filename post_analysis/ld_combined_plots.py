"""
Combined 3-panel "mountain" plot per hit, all panels sharing the SNP-index x-axis:

  panel 1 (top)    : global (all-samples) LD r^2 triangle, CLIPPED at r^2 = 0.25
                     (vmax=0.25), with vertical guides at EVERY window boundary
                     (significant windows + ±6 extension windows).
  panel 2 (middle) : per-locus haplotype entropy (±7 SNP window, all samples).
  panel 3 (bottom) : association track — bars per SNP, grey with the fine-mapped
                     concordant SNPs coloured (red risk / blue protective).

Two versions per hit (identical top two panels, only the bottom track changes):
  {hit}__combined_rawP.png  -> bottom = -log10(raw p)
  {hit}__combined_BYp.png   -> bottom = -log10(BY-corrected p)

Heavy: extracts phased haplotypes from the UKBB VCF (bcftools slice + plink2
--export haps) to compute the entropy, exactly like ld_entropy_plots.py.

Output: ld_mountain_plots_all/{hit}__combined_rawP.* and __combined_BYp.*
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
from matplotlib.patches import Patch
from statsmodels.stats.multitest import multipletests

BASE     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb'
LD_DIR   = f'{BASE}/SuSiE_ld_6wind'
FM_DIRS  = [f'{BASE}/fine_mapping_new', f'{BASE}/fine_mapping_new_6wind']
CAND_FILE = f'{BASE}/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
VCF      = f'{BASE}/vcf_file/ukbb.vcf.gz'
KCUTOFF  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
PLINK2   = '/private/home/rsmerigl/plink2'
WIND_COVAR_DIRS = ['/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new',
                   '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new_6wind']
OUT_DIR  = os.path.join(os.path.dirname(__file__), 'ld_mountain_plots_aggregate')
WORK     = os.path.join(OUT_DIR, '_entropy_work')
# shared, persistent entropy cache (computed once, reused by every script)
ENTROPY_CACHE = f'{BASE}/entropy_cache'

# Crop the LD triangle vertically to a TRAPEZOID: keep only the near-diagonal
# band (short/medium-range LD, where the signal is) and drop the empty apex.
# Fraction of the full triangle height (N/2) that stays visible. Tune freely.
TRIANGLE_HEIGHT_FRAC = 0.40
HALF    = 7          # ±7 SNP entropy window (15 SNPs)
C_OTHER, C_RISK, C_PROT = '#bdbdbd', '#d62728', '#1f77b4'

mpl.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans']})

_cand = pd.read_csv(CAND_FILE, sep='\t')
_cand = _cand[_cand['concordant_direction'] == True]                     # noqa: E712
CONCORDANT_POS = {h: set(g['POS'].astype(int)) for h, g in _cand.groupby('hit')}


# ── fine-mapping stats ─────────────────────────────────────────────────────────
def load_fm(hit):
    rows = []
    for fm_dir in FM_DIRS:
        for f in glob.glob(os.path.join(fm_dir, hit, '*.glm.logistic.hybrid')):
            d = pd.read_csv(f, sep='\t', dtype=str, low_memory=False)
            rows.append(d[d['TEST'] == 'ADD'][['ID', 'POS', 'P', 'OR']])
    if not rows:
        return None
    d = pd.concat(rows, ignore_index=True).drop_duplicates('ID')
    for c in ('POS', 'P', 'OR'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    return d.set_index('ID')


def window_boundaries(hit):
    bnds = set()
    prefix = f'{hit}_'
    for wdir in WIND_COVAR_DIRS:
        if not os.path.isdir(wdir):
            continue
        for fn in os.listdir(wdir):
            if fn.startswith(prefix) and fn.endswith('_covar.tsv'):
                parts = fn.replace('_covar.tsv', '').split('_')
                if len(parts) >= 5:
                    bnds.add(int(parts[3])); bnds.add(int(parts[4]))
    return sorted(bnds)


def boundary_x(pos, boundaries):
    m = np.isfinite(pos)
    if m.sum() < 2:
        return []
    order = np.argsort(pos[m])
    xp = pos[m][order]; fp = np.arange(len(pos))[m][order]
    lo, hi = xp.min(), xp.max()
    return [float(np.interp(b, xp, fp)) for b in boundaries if lo <= b <= hi]


# ── haplotype extraction + entropy (all samples), with a persistent cache ──────
def cached_entropy(hit, snps, pos):
    """Return per-SNP entropy for a hit, reading the shared cache if present
       (keyed by SNP ID) and otherwise computing it once and saving it."""
    cf = os.path.join(ENTROPY_CACHE, f'{hit}_entropy.tsv')
    if os.path.exists(cf):
        d = pd.read_csv(cf, sep='\t').set_index('ID')['entropy']
        print(f'    entropy: reusing cache {os.path.basename(cf)}')
        return d.reindex(snps).to_numpy(dtype=float)
    ent = compute_entropy(hit, snps, pos)
    os.makedirs(ENTROPY_CACHE, exist_ok=True)
    pd.DataFrame({'ID': snps, 'entropy': ent}).to_csv(cf, sep='\t', index=False)
    print(f'    entropy: computed and cached -> {os.path.basename(cf)}')
    return ent


def compute_entropy(hit, snps, pos):
    finite = pos[np.isfinite(pos)]
    if finite.size == 0:
        return np.full(len(snps), np.nan)
    chrom = hit.split('_')[-1].replace('chr', '')
    region = f'{chrom}:{int(finite.min())}-{int(finite.max())}'
    os.makedirs(WORK, exist_ok=True)
    region_vcf = os.path.join(WORK, f'{hit}_region.vcf.gz')
    out_base   = os.path.join(WORK, hit)
    var_file   = os.path.join(LD_DIR, hit, f'{hit}_ld_6wind.phased.vcor1.vars')
    try:
        subprocess.run(['bcftools', 'view', '-r', region, VCF, '-Oz', '-o', region_vcf], check=True)
        subprocess.run(['tabix', '-f', '-p', 'vcf', region_vcf], check=True)
        subprocess.run([PLINK2, '--vcf', region_vcf, '--keep', KCUTOFF,
                        '--extract', var_file, '--export', 'haps', '--out', out_base],
                       check=True, stdout=subprocess.DEVNULL)
        ids, rows = [], []
        with open(out_base + '.haps') as f:
            for line in f:
                p = line.split(' ', 5)
                ids.append(p[1])
                rows.append((np.frombuffer(p[5].rstrip('\n').encode(), np.uint8)[::2] - 48).astype(np.uint8))
        H = np.vstack(rows)
        idpos = {v: i for i, v in enumerate(ids)}
        order = [idpos.get(v, -1) for v in snps]
        present = np.array([o >= 0 for o in order])
        Hal = np.vstack([H[o] if o >= 0 else np.zeros(H.shape[1], np.uint8) for o in order])
        V, twoN = Hal.shape
        ent = np.full(V, np.nan)
        for i in range(V):
            idx = [j for j in range(max(0, i - HALF), min(V, i + HALF + 1)) if present[j]]
            if len(idx) < 2:
                continue
            codes = (1 << np.arange(len(idx))).astype(np.int64) @ Hal[idx]
            _, c = np.unique(codes, return_counts=True)
            pf = c / twoN
            ent[i] = -(pf * np.log2(pf)).sum()
        return ent
    finally:
        for pth in (region_vcf, region_vcf + '.tbi', out_base + '.haps',
                    out_base + '.sample', out_base + '.log'):
            if os.path.exists(pth):
                os.remove(pth)


# ── figure ─────────────────────────────────────────────────────────────────────
def draw_panels(hit, track, ld2, N, pos, P, ORv, ent, bxs):
    # triangle panel height scales with the crop fraction so the trapezoid keeps
    # the same cell proportions as the full pyramid (not vertically elongated)
    tri_h = 3.0 * TRIANGLE_HEIGHT_FRAC
    total = tri_h + 1.2 + 1.4 + 1.4
    fig = plt.figure(figsize=(9, 1.3 * total))
    gs = fig.add_gridspec(4, 2, height_ratios=[tri_h, 1.2, 1.4, 1.4], width_ratios=[45, 1],
                          hspace=0.07, wspace=0.02)
    axm = fig.add_subplot(gs[0, 0])
    axe = fig.add_subplot(gs[1, 0], sharex=axm)
    axp = fig.add_subplot(gs[2, 0], sharex=axm)
    axo = fig.add_subplot(gs[3, 0], sharex=axm)
    cax = fig.add_subplot(gs[0, 1])

    # panel 1: LD triangle clipped at r^2 = 0.25
    ii, jj = np.meshgrid(np.arange(N + 1), np.arange(N + 1), indexing='ij')
    X = (ii + jj) / 2.0; Y = (jj - ii) / 2.0
    C = ld2.copy().astype(float); C[np.tril_indices(N, -1)] = np.nan
    pc = axm.pcolormesh(X, Y, C, cmap='Reds', vmin=0, vmax=1, shading='flat', rasterized=True)
    axm.set_ylim(0, TRIANGLE_HEIGHT_FRAC * N / 2.0)     # horizontal crop -> trapezoid
    axm.set_xlim(0, N); axm.set_yticks([])
    axm.tick_params(axis='x', labelbottom=False, bottom=False)
    for s in ('top', 'right', 'left'):
        axm.spines[s].set_visible(False)
    axm.set_title(f"{hit}  (all samples LD, ±6 windows)",
                  fontsize=9, fontweight='bold')
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r'LD $r^2$', fontsize=7); cb.ax.tick_params(labelsize=6)

    # panel 2: haplotype entropy
    x = np.arange(N)
    axe.fill_between(x, 0, np.nan_to_num(ent, nan=0.0), where=np.isfinite(ent),
                     color='#6a3d9a', alpha=0.25, linewidth=0)
    axe.plot(x, ent, color='#6a3d9a', lw=0.6, rasterized=True)
    axe.set_ylabel('Haplotype entropy\n(±7 SNP, bits)', fontsize=8)
    axe.set_ylim(bottom=0)
    axe.tick_params(axis='x', labelbottom=False, bottom=False)
    axe.tick_params(axis='y', labelsize=6)
    for s in ('top', 'right'):
        axe.spines[s].set_visible(False)

    # panel 3: association track (raw p or BY p)
    if track == 'BYp':
        finite = np.isfinite(P) & (P > 0)
        adj = np.full(N, np.nan)
        if finite.sum() > 0:
            _, p_adj, _, _ = multipletests(P[finite], alpha=0.05, method='fdr_by')
            adj[finite] = p_adj
        height = -np.log10(np.where((adj > 0) & np.isfinite(adj), adj, np.nan))
        ylab = r'$-\log_{10}(\mathrm{BY}\ p)$'
    else:
        height = -np.log10(np.where((P > 0) & np.isfinite(P), P, np.nan))
        ylab = r'$-\log_{10}(p)$'

    cand_pos = CONCORDANT_POS.get(hit, set())
    is_cand = np.array([(np.isfinite(p) and int(p) in cand_pos) for p in pos])
    colors = np.full(N, C_OTHER, dtype=object)
    colors[is_cand & (ORv >= 1)] = C_RISK
    colors[is_cand & (ORv <  1)] = C_PROT
    valid = np.isfinite(height)
    other = valid & ~is_cand
    axp.bar(x[other], height[other], width=1.0, color=C_OTHER, linewidth=0, rasterized=True)
    hicol = valid & is_cand
    axp.bar(x[hicol], height[hicol], width=1.0, color=list(colors[hicol]),
            linewidth=0, rasterized=True, zorder=3)
    axp.set_ylabel(ylab, fontsize=8); axp.set_ylim(bottom=0)
    axp.tick_params(axis='y', labelsize=6)
    axp.tick_params(axis='x', labelbottom=False, bottom=False)   # OR panel is the bottom one now
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)
    handles = [Patch(facecolor=C_OTHER, label='Other SNPs')]
    if (hicol & (ORv >= 1)).any():
        handles.append(Patch(facecolor=C_RISK, label='Fine-mapped · concordant (risk)'))
    if (hicol & (ORv < 1)).any():
        handles.append(Patch(facecolor=C_PROT, label='Fine-mapped · concordant (protective)'))
    axp.legend(handles=handles, fontsize=6, loc='upper right', frameon=False,
               handlelength=1.0, borderaxespad=0.2)

    # panel 4: signed OR bars around the OR = 1 line (red = risk, blue = protective)
    valid_or = np.isfinite(ORv) & (ORv > 0)
    or_colors = np.where(ORv[valid_or] >= 1, C_RISK, C_PROT)
    axo.bar(x[valid_or], ORv[valid_or] - 1, bottom=1.0, width=1.0, color=or_colors,
            linewidth=0, rasterized=True)
    axo.axhline(1.0, color='black', lw=0.5)
    axo.set_ylabel('OR', fontsize=8)
    if valid_or.any():
        lo, hi = np.nanpercentile(ORv[valid_or], [1, 99])
        axo.set_ylim(min(0.9, lo), max(1.1, hi))
    axo.tick_params(axis='y', labelsize=6)
    for s in ('top', 'right'):
        axo.spines[s].set_visible(False)
    or_handles = [Patch(facecolor=C_RISK, label='OR ≥ 1 (risk)'),
                  Patch(facecolor=C_PROT, label='OR < 1 (protective)')]
    axo.legend(handles=or_handles, fontsize=6, loc='upper right', frameon=False,
               handlelength=1.0, borderaxespad=0.2)

    # vertical window guides ONLY on the pyramid (triangle) panel
    for xb in bxs:
        axm.axvline(xb, color='#2b2b2b', lw=0.35, ls=(0, (3, 3)), alpha=0.45, zorder=4)

    nt = 6
    idx = np.linspace(0, N - 1, nt).astype(int)
    labs = [f'{pos[i]/1e6:.2f}' if np.isfinite(pos[i]) else '' for i in idx]
    axo.set_xticks(idx); axo.set_xticklabels(labs, fontsize=6)
    chrom = hit.split('_')[-1]
    axo.set_xlabel(f'Position on {chrom} (Mb)', fontsize=8)

    os.makedirs(OUT_DIR, exist_ok=True)
    suffix = 'BYp' if track == 'BYp' else 'rawP'
    base = os.path.join(OUT_DIR, f'{hit}__combined_{suffix}')
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{base}.pdf', dpi=300, bbox_inches='tight', format='pdf')
    plt.close(fig)
    print(f'    saved -> {os.path.basename(base)}.png/.pdf')


def make_plots(hit):
    ld_file  = os.path.join(LD_DIR, hit, f'{hit}_ld_6wind.phased.vcor1')
    var_file = ld_file + '.vars'
    if not os.path.exists(ld_file):
        print(f'  {hit}: no LD matrix, skip'); return
    snps = pd.read_csv(var_file, header=None)[0].astype(str).tolist()
    N = len(snps)
    fm = load_fm(hit)
    if fm is None:
        print(f'  {hit}: no fine-mapping stats, skip'); return
    pos = fm['POS'].reindex(snps).to_numpy(dtype=float)
    P   = fm['P'].reindex(snps).to_numpy(dtype=float)
    ORv = fm['OR'].reindex(snps).to_numpy(dtype=float)

    print(f'  {hit}: {N} SNPs — entropy ...')
    ent = cached_entropy(hit, snps, pos)
    R = np.loadtxt(ld_file); ld2 = R ** 2
    bxs = boundary_x(pos, window_boundaries(hit))

    for track in ('rawP', 'BYp'):
        draw_panels(hit, track, ld2, N, pos, P, ORv, ent, bxs)


if __name__ == '__main__':
    import sys
    hits = sys.argv[1:] or sorted(
        d for d in os.listdir(LD_DIR) if os.path.isdir(os.path.join(LD_DIR, d)))
    for hit in hits:
        try:
            make_plots(hit)
        except Exception as e:
            print(f'  {hit}: FAILED ({e})')
    if os.path.isdir(WORK) and not os.listdir(WORK):
        shutil.rmtree(WORK, ignore_errors=True)
