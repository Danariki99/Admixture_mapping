"""
Unified "mountain" / TAD-style plots per hit (EXTENDED, ±6 windows).

Produces, for every hit, the 6 LD-variant plots, all in a SINGLE folder with a
coherent naming scheme  {hit}__{ldsource}__{track}.png/.pdf :

  ldsource:                          track (bottom panel):
    allsamples  = SuSiE_ld_6wind        OR   = signed OR bars around OR=1 (red risk / blue prot.)
                  (KING-cutoff cohort,  BYp  = -log10(BY-corrected p), grey + fine-mapped-concordant coloured
                   all ancestries)      P    = -log10(raw p),          grey + fine-mapped-concordant coloured
    hgdp        = SuSiE_ld_HGDP_6wind
                  (HGDP reference panel of the hit's ancestry = "non-admixed")

The top panel is always the LD r^2 triangle for the chosen LD source; the bottom
panel is the association track (from fine_mapping_new + fine_mapping_new_6wind).

NOTE on SNP matching: the all-samples LD .vars carry rsIDs, the HGDP .vars carry
positional IDs (chr:pos). Everything is therefore matched to the fine-mapping
stats and to the concordant candidates BY POSITION.

The entropy plot ({hit}__allsamples__entropy) is produced by the companion
script ld_entropy_plots.py (heavy, needs haplotype extraction).

Output: ld_mountain_plots_all/{hit}__{ldsource}__{track}.png / .pdf
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from statsmodels.stats.multitest import multipletests

BASE     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb'
FM_DIRS  = [f'{BASE}/fine_mapping_new', f'{BASE}/fine_mapping_new_6wind']
CAND_FILE = f'{BASE}/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
# window covar folders: significant windows + the ±6 extension windows.
# ALL of them must be delimited (not only the 6+6 extensions).
WIND_COVAR_DIRS = ['/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new',
                   '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new_6wind']
OUT_DIR  = os.path.join(os.path.dirname(__file__), 'ld_mountain_plots_all')

# ld source -> (folder, ld-file suffix, id-type of the .vars entries)
LD_SOURCES = {
    'allsamples': (f'{BASE}/SuSiE_ld_6wind',      '_ld_6wind',      'rsid', 'all samples'),
    'hgdp':       (f'{BASE}/SuSiE_ld_HGDP_6wind', '_ld_HGDP_6wind', 'pos',  'HGDP'),
}
TRACKS = ['OR', 'BYp', 'P']

C_OTHER, C_RISK, C_PROT = '#bdbdbd', '#d62728', '#1f77b4'

mpl.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans']})

# ── fine-mapping candidates concordant with the admixture signal (POS + OR) ─────
_cand = pd.read_csv(CAND_FILE, sep='\t')
_cand = _cand[_cand['concordant_direction'] == True]                     # noqa: E712
CONCORDANT_POS = {h: set(g['POS'].astype(int)) for h, g in _cand.groupby('hit')}


def load_fm(hit):
    """Fine-mapping ADD stats for a hit. Returns two lookups:
       by_id : DataFrame indexed by rsID  with columns POS, P, OR
       by_pos: DataFrame indexed by POS    with columns P, OR (one row per POS)."""
    rows = []
    for fm_dir in FM_DIRS:
        for f in glob.glob(os.path.join(fm_dir, hit, '*.glm.logistic.hybrid')):
            d = pd.read_csv(f, sep='\t', dtype=str, low_memory=False)
            rows.append(d[d['TEST'] == 'ADD'][['ID', 'POS', 'P', 'OR']])
    if not rows:
        return None, None
    d = pd.concat(rows, ignore_index=True).drop_duplicates('ID')
    for c in ('POS', 'P', 'OR'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    by_id  = d.set_index('ID')
    by_pos = d.dropna(subset=['POS']).drop_duplicates('POS').set_index('POS')[['P', 'OR']]
    return by_id, by_pos


def window_boundaries(hit):
    """Sorted unique start/end positions of ALL windows of a hit
       (significant windows + ±6 extension windows)."""
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
    """Map genomic boundary positions to fractional SNP-index x-coordinates."""
    m = np.isfinite(pos)
    if m.sum() < 2:
        return []
    order = np.argsort(pos[m])
    xp = pos[m][order]
    fp = np.arange(len(pos))[m][order]
    lo, hi = xp.min(), xp.max()
    return [float(np.interp(b, xp, fp)) for b in boundaries if lo <= b <= hi]


def snp_pos_stats(snps, idtype, by_id, by_pos):
    """Return (pos, P, OR) arrays aligned to snps, matched appropriately per source."""
    if idtype == 'rsid':
        pos = by_id['POS'].reindex(snps).to_numpy(dtype=float)
        P   = by_id['P'].reindex(snps).to_numpy(dtype=float)
        OR  = by_id['OR'].reindex(snps).to_numpy(dtype=float)
    else:  # positional ids "chr:pos" -> match fine mapping by POS
        pos = np.array([int(s.split(':')[1]) if ':' in s else -1 for s in snps], dtype=float)
        P   = by_pos['P'].reindex(pos).to_numpy(dtype=float)
        OR  = by_pos['OR'].reindex(pos).to_numpy(dtype=float)
    return pos, P, OR


def draw_triangle(fig, axm, cax, ld2, N, title):
    ii, jj = np.meshgrid(np.arange(N + 1), np.arange(N + 1), indexing='ij')
    X = (ii + jj) / 2.0
    Y = (jj - ii) / 2.0
    C = ld2.copy().astype(float)
    C[np.tril_indices(N, -1)] = np.nan          # upper triangle only
    pc = axm.pcolormesh(X, Y, C, cmap='Reds', vmin=0, vmax=1, shading='flat', rasterized=True)
    axm.set_ylim(0, N / 2.0)
    axm.set_xlim(0, N)
    axm.set_yticks([])
    axm.tick_params(axis='x', labelbottom=False, bottom=False)
    for s in ('top', 'right', 'left'):
        axm.spines[s].set_visible(False)
    axm.set_title(title, fontsize=9, fontweight='bold')
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r'LD $r^2$', fontsize=7); cb.ax.tick_params(labelsize=6)


def make_plot(hit, ldsource, track, by_id, by_pos):
    folder, suffix, idtype, srclabel = LD_SOURCES[ldsource]
    ld_file  = os.path.join(folder, hit, f'{hit}{suffix}.phased.vcor1')
    var_file = ld_file + '.vars'
    if not os.path.exists(ld_file):
        print(f'  {hit} [{ldsource}/{track}]: no LD matrix, skip'); return

    snps = pd.read_csv(var_file, header=None)[0].astype(str).tolist()
    N = len(snps)
    R = np.loadtxt(ld_file)
    ld2 = R ** 2

    pos, P, ORv = snp_pos_stats(snps, idtype, by_id, by_pos)

    fig = plt.figure(figsize=(9, 5.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.4], width_ratios=[45, 1],
                          hspace=0.05, wspace=0.02)
    axm = fig.add_subplot(gs[0, 0])
    axp = fig.add_subplot(gs[1, 0], sharex=axm)
    cax = fig.add_subplot(gs[0, 1])

    title = f"{hit.replace('NAT', 'AMR')}  ({srclabel} LD, ±6 windows)"
    draw_triangle(fig, axm, cax, ld2, N, title)

    x = np.arange(N)
    cand_pos = CONCORDANT_POS.get(hit, set())
    is_cand = np.array([(np.isfinite(p) and int(p) in cand_pos) for p in pos])

    if track == 'OR':
        valid = np.isfinite(ORv) & (ORv > 0)
        colors = np.where(ORv[valid] >= 1, C_RISK, C_PROT)
        axp.bar(x[valid], ORv[valid] - 1, bottom=1.0, width=1.0, color=colors,
                linewidth=0, rasterized=True)
        axp.axhline(1.0, color='black', lw=0.5)
        axp.set_ylabel('OR', fontsize=8)
        if valid.any():
            lo, hi = np.nanpercentile(ORv[valid], [1, 99])
            axp.set_ylim(min(0.9, lo), max(1.1, hi))
        handles = [Patch(facecolor=C_RISK, label='OR ≥ 1 (risk)'),
                   Patch(facecolor=C_PROT, label='OR < 1 (protective)')]
        axp.legend(handles=handles, fontsize=6, loc='upper right', frameon=False,
                   handlelength=1.0, borderaxespad=0.2)
    else:
        if track == 'BYp':
            finite = np.isfinite(P) & (P > 0)
            adj = np.full(N, np.nan)
            if finite.sum() > 0:
                _, p_adj, _, _ = multipletests(P[finite], alpha=0.05, method='fdr_by')
                adj[finite] = p_adj
            height = -np.log10(np.where((adj > 0) & np.isfinite(adj), adj, np.nan))
            ylab = r'$-\log_{10}(\mathrm{BY}\ p)$'
        else:  # raw P
            height = -np.log10(np.where((P > 0) & np.isfinite(P), P, np.nan))
            ylab = r'$-\log_{10}(p)$'

        colors = np.full(N, C_OTHER, dtype=object)
        colors[is_cand & (ORv >= 1)] = C_RISK
        colors[is_cand & (ORv <  1)] = C_PROT
        valid = np.isfinite(height)
        other = valid & ~is_cand
        axp.bar(x[other], height[other], width=1.0, color=C_OTHER, linewidth=0, rasterized=True)
        hicol = valid & is_cand
        axp.bar(x[hicol], height[hicol], width=1.0, color=list(colors[hicol]),
                linewidth=0, rasterized=True, zorder=3)
        axp.set_ylabel(ylab, fontsize=8)
        axp.set_ylim(bottom=0)
        handles = [Patch(facecolor=C_OTHER, label='Other SNPs')]
        if (hicol & (ORv >= 1)).any():
            handles.append(Patch(facecolor=C_RISK, label='Fine-mapped · concordant (risk)'))
        if (hicol & (ORv < 1)).any():
            handles.append(Patch(facecolor=C_PROT, label='Fine-mapped · concordant (protective)'))
        axp.legend(handles=handles, fontsize=6, loc='upper right', frameon=False,
                   handlelength=1.0, borderaxespad=0.2)

    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)

    # vertical guides at the window boundaries — ONLY on the pyramid (triangle)
    for xb in boundary_x(pos, window_boundaries(hit)):
        axm.axvline(xb, color='#2b2b2b', lw=0.35, ls=(0, (3, 3)), alpha=0.45, zorder=4)

    nt = 6
    idx = np.linspace(0, N - 1, nt).astype(int)
    labs = [f'{pos[i]/1e6:.2f}' if np.isfinite(pos[i]) else '' for i in idx]
    axp.set_xticks(idx); axp.set_xticklabels(labs, fontsize=6)
    chrom = hit.split('_')[-1]
    axp.set_xlabel(f'Position on {chrom} (Mb)', fontsize=8)
    axp.tick_params(axis='y', labelsize=6)

    os.makedirs(OUT_DIR, exist_ok=True)
    base = os.path.join(OUT_DIR, f'{hit}__{ldsource}__{track}')
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{base}.pdf', dpi=300, bbox_inches='tight', format='pdf')
    plt.close(fig)
    print(f'    saved -> {os.path.basename(base)}.png/.pdf')


if __name__ == '__main__':
    import sys
    hits = sys.argv[1:] or sorted(
        d for d in os.listdir(LD_SOURCES['allsamples'][0])
        if os.path.isdir(os.path.join(LD_SOURCES['allsamples'][0], d)))
    for hit in hits:
        print(f'{hit}: loading fine-mapping stats ...')
        by_id, by_pos = load_fm(hit)
        if by_id is None:
            print('  no fine-mapping stats, skip'); continue
        for ldsource in LD_SOURCES:
            for track in TRACKS:
                make_plot(hit, ldsource, track, by_id, by_pos)
