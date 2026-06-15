"""
Two figures from the annotated fine-mapping candidate SNPs:

  1. fine_mapping_OR_barplot.{png,pdf}
     Grouped bar plot. Y = OR (bars emanate from the OR=1 null line). X = hits
     (ancestry x phenotype), each a group; within a group one bar per SNP.
     Bars coloured by functional consequence category (Intronic, 3' UTR, ...).

  2. fine_mapping_consequence_pie.{png,pdf}
     Pie chart of the consequence-category composition, SAME palette.

OR is hit-specific (comes from the per-hit fine-mapping regression), so the bar
plot uses the full 179 (SNP x hit) rows from fine_mapping_all_candidates.tsv.
The consequence category is a property of the variant, mapped by SNP ID from the
annotated file. The pie chart uses the unique SNPs (one per ID).

Inputs:
  fine_mapping_all_candidates.tsv        (cols: hit, ID, OR, beta, ...)  -> per-hit OR
  fine_mapping_candidates_annotated.tsv  (cols: ID, category, ...)       -> consequence
"""

import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

CAND_FILE  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
ANNOT_FILE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_candidates_annotated.tsv'
EXCEL_PATH = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx'
OUT_DIR    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/plots/ukbb'

MM_TO_INCH = 1 / 25.4

mpl.rcParams.update({
    'font.size':       6,
    'axes.titlesize':  7,
    'axes.labelsize':  6,
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'pdf.fonttype':    42,
    'ps.fonttype':     42,
    'savefig.dpi':     600,
    'font.family':     'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
})

# Fixed palette + severity order so both figures match and colours are stable.
CATEGORY_ORDER = [
    'Stop gained', 'Missense', 'Synonymous', 'Non-coding exon',
    "5' UTR", "3' UTR", 'Intronic', 'Upstream', 'Downstream', 'Intergenic',
]
CATEGORY_COLORS = {
    'Stop gained':     '#d62728',
    'Missense':        '#ff7f0e',
    'Synonymous':      '#bcbd22',
    'Non-coding exon': '#8c564b',
    "5' UTR":          '#1f77b4',
    "3' UTR":          '#17becf',
    'Intronic':        '#2ca02c',
    'Upstream':        '#9467bd',
    'Downstream':      '#e377c2',
    'Intergenic':      '#7f7f7f',
}
FALLBACK_COLOR = '#cccccc'


def color_for(cat):
    return CATEGORY_COLORS.get(cat, FALLBACK_COLOR)


def load_excel_labels(path):
    try:
        return pd.read_excel(path, sheet_name='first_batch', usecols='B:C')
    except Exception as e:
        print(f'  (excel labels unavailable: {e})')
        return None


def pheno_label(excel_df, pheno):
    if excel_df is None:
        label = pheno
    else:
        row = excel_df.loc[excel_df['ID'] == pheno, 'ID2']
        label = row.iloc[0] if not row.empty else pheno
    # drop the leading code (split on '_', remove the first token)
    parts = label.split('_')
    if len(parts) > 1:
        label = '_'.join(parts[1:])
    return label.replace('_', ' ')


def hit_label(hit, excel_df):
    # AFR_HC219_chr6 -> "AFR – <phenotype name>"  (single line, for angled labels)
    parts = hit.split('_')
    anc   = parts[0]
    pheno = parts[1] if len(parts) > 1 else ''
    name  = pheno_label(excel_df, pheno)
    return f'{anc} – {name}' if name else anc


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    excel_df = load_excel_labels(EXCEL_PATH)

    # per-hit OR (179 SNP x hit rows)
    df = pd.read_csv(CAND_FILE, sep='\t')
    df['OR'] = pd.to_numeric(df['OR'], errors='coerce')
    df = df.dropna(subset=['OR']).copy()

    # consequence category mapped by SNP ID (variant property)
    annot = pd.read_csv(ANNOT_FILE, sep='\t')
    id2cat = dict(zip(annot['ID'], annot['category']))
    df['category'] = df['ID'].map(id2cat).fillna('Unannotated')
    print(f'Loaded {len(df)} SNP x hit rows across {df["hit"].nunique()} hits '
          f'({df["ID"].nunique()} unique SNPs)')

    # categories actually present, in severity order
    present = [c for c in CATEGORY_ORDER if c in set(df['category'])]
    extra   = [c for c in df['category'].dropna().unique() if c not in CATEGORY_ORDER]
    present = present + sorted(extra)

    # ── Figure 1: grouped bar plot ────────────────────────────────────────────
    # order: asthma hits first, then hypothyroidism (shorter labels end up on the
    # right, saving space near the legend)
    def _order(hit):
        p  = hit.split('_')
        nm = pheno_label(excel_df, p[1] if len(p) > 1 else '').lower()
        return (0 if 'asthma' in nm else 1, hit)
    hits = sorted(df['hit'].unique(), key=_order)
    GAP  = 6.0          # gap (in bar units) between hit groups
    x = 0.0
    bar_x, bar_h, bar_c = [], [], []
    group_centers, group_labels = [], []
    group_lead = []   # (x, OR, rsID) of the strongest-effect SNP per hit

    for hit in hits:
        sub = df[df['hit'] == hit].sort_values('OR', ascending=False)
        start = x
        g_x, g_or, g_id = [], [], []
        for _, row in sub.iterrows():
            bar_x.append(x)
            bar_h.append(row['OR'])
            bar_c.append(color_for(row['category']))
            g_x.append(x); g_or.append(row['OR']); g_id.append(row['ID'])
            x += 1.0
        end = x
        # lead SNP = strongest effect = max |log(OR)| (furthest from the null)
        lead_i = int(np.argmax(np.abs(np.log(np.array(g_or)))))
        group_lead.append((g_x[lead_i], g_or[lead_i], g_id[lead_i]))
        group_centers.append((start + end - 1) / 2)
        group_labels.append(hit_label(hit, excel_df))
        x += GAP

    fig, ax = plt.subplots(figsize=(183 * MM_TO_INCH, 75 * MM_TO_INCH))
    # bars emanate from the OR = 1 null line
    ax.bar(bar_x, [h - 1 for h in bar_h], bottom=1.0, width=0.9,
           color=bar_c, linewidth=0)
    ax.axhline(1.0, color='black', linewidth=0.5, zorder=3)

    # rsID of the lead SNP: above its bar tilted right (risk) /
    #                       below its bar tilted left (protective)
    for gx, gor, gid in group_lead:
        if gor >= 1:
            ax.annotate(gid, (gx, gor), xytext=(0, 2), textcoords='offset points',
                        rotation=45, ha='left', va='bottom', fontsize=4)
        else:
            ax.annotate(gid, (gx, gor), xytext=(0, -2), textcoords='offset points',
                        rotation=45, ha='right', va='top', fontsize=4)

    ax.set_xticks(group_centers)
    ax.set_xticklabels(group_labels, fontsize=5, rotation=35, ha='right')
    ax.set_xlim(-1, x - GAP)
    ax.margins(y=0.12)   # headroom for the rsID labels
    ax.set_ylabel('Odds ratio (OR)', fontsize=6)
    ax.set_title('Fine-mapping candidate SNPs', fontsize=7)
    ax.tick_params(axis='both', width=0.5, length=2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

    # legend on the right, one category per line (fully readable names)
    handles = [Patch(facecolor=color_for(c), label=c) for c in present]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.01, 0.5),
              ncol=1, frameon=False, fontsize=5, handletextpad=0.5,
              labelspacing=0.4, title='Consequence', title_fontsize=6)

    plt.tight_layout()
    base1 = os.path.join(OUT_DIR, 'fine_mapping_OR_barplot')
    plt.savefig(f'{base1}.png', bbox_inches='tight', pad_inches=0.02)
    plt.savefig(f'{base1}.pdf', bbox_inches='tight', pad_inches=0.02, format='pdf')
    plt.close()
    print(f'Saved → {base1}.png / .pdf')

    # ── Figure 2: consequence pie chart (all 179 SNP x hit rows) ───────────────
    counts = df['category'].value_counts()
    counts = counts.reindex([c for c in present if c in counts.index])
    pie_colors = [color_for(c) for c in counts.index]

    fig, ax = plt.subplots(figsize=(90 * MM_TO_INCH, 75 * MM_TO_INCH))
    wedges, _texts, autotexts = ax.pie(
        counts.values, colors=pie_colors, startangle=90, counterclock=False,
        autopct=lambda p: f'{p:.0f}%' if p >= 4 else '',
        pctdistance=0.75, wedgeprops={'linewidth': 0.3, 'edgecolor': 'white'},
    )
    for t in autotexts:
        t.set_fontsize(5)
    ax.set_title('Consequence composition of candidate SNPs', fontsize=7)
    ax.axis('equal')

    handles = [Patch(facecolor=color_for(c), label=f'{c} ({counts[c]})')
               for c in counts.index]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
              frameon=False, fontsize=5)

    plt.tight_layout()
    base2 = os.path.join(OUT_DIR, 'fine_mapping_consequence_pie')
    plt.savefig(f'{base2}.png', bbox_inches='tight', pad_inches=0.02)
    plt.savefig(f'{base2}.pdf', bbox_inches='tight', pad_inches=0.02, format='pdf')
    plt.close()
    print(f'Saved → {base2}.png / .pdf')


if __name__ == '__main__':
    main()
