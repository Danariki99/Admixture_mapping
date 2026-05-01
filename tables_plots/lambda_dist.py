import os
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

HITS = [
    ('EUR', 'HC1007'), ('EUR', 'HC1581'), ('EUR', 'HC221'), ('EUR', 'HC382'),
    ('SAS', 'FH1220'), ('SAS', 'HC1036'), ('SAS', 'HC1581'), ('SAS', 'HC382'),
    ('WAS', 'HC219'),  ('WAS', 'HC221'),  ('WAS', 'HC643'),
]
ANCESTRIES = ['EUR', 'SAS', 'WAS']

def glob_logs(base_folder):
    for root, dirs, files in os.walk(base_folder):
        for f in files:
            if f == 'output.log':
                yield os.path.join(root, f)

def parse_lambda(log_path):
    with open(log_path) as f:
        for line in f:
            m = re.search(r'lambda \(based on median chisq\) = ([\d.]+?)\.?\s', line)
            if m:
                return float(m.group(1))
    return None

def get_anc_pheno(log_path):
    m = re.search(r'output_ancestry_(\w+)/(\w+)/', log_path)
    if m:
        return m.group(1), m.group(2)
    return None, None

def extract_lambdas(base_folder, hits_only=False, ancestry=None):
    lambdas = []
    for log in glob_logs(base_folder):
        anc, pheno = get_anc_pheno(log)
        if hits_only and (anc, pheno) not in HITS:
            continue
        if ancestry and anc != ancestry:
            continue
        val = parse_lambda(log)
        if val is not None:
            lambdas.append(val)
    return lambdas

def violin(ax, data, pos, color):
    p = ax.violinplot([data], positions=[pos], showmedians=True, showextrema=True)
    p['bodies'][0].set_facecolor(color)
    p['bodies'][0].set_alpha(0.7)
    p['cmedians'].set_color('black')
    p['cmedians'].set_linewidth(2)


ALL_ANCESTRIES = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']

old_folder    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output_old/ukbb'
cutoff_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output_distributions_kcutoff_177/ukbb'
new_folder    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb'

# --- Data ---
lambdas_old_all  = extract_lambdas(old_folder,    hits_only=False)
lambdas_old_hits = extract_lambdas(old_folder,    hits_only=True)
lambdas_cutoff   = extract_lambdas(cutoff_folder, hits_only=True)
lambdas_new_all  = extract_lambdas(new_folder,    hits_only=False)

lambdas_by_anc_hits = {
    anc: {
        'old':    extract_lambdas(old_folder,    hits_only=True, ancestry=anc),
        'cutoff': extract_lambdas(cutoff_folder, hits_only=True, ancestry=anc),
    }
    for anc in ANCESTRIES
}

lambdas_by_anc_all = {
    anc: {
        'old': extract_lambdas(old_folder,  hits_only=False, ancestry=anc),
        'new': extract_lambdas(new_folder,  hits_only=False, ancestry=anc),
    }
    for anc in ALL_ANCESTRIES
}

print(f"Original complete:  n={len(lambdas_old_all)},  median={np.median(lambdas_old_all):.3f}")
print(f"Original hits:      n={len(lambdas_old_hits)}, median={np.median(lambdas_old_hits):.3f}")
print(f"Cutoff 177 hits:    n={len(lambdas_cutoff)},   median={np.median(lambdas_cutoff):.3f}")
print(f"New analysis (all): n={len(lambdas_new_all)},  median={np.median(lambdas_new_all):.3f}")

color_old    = 'steelblue'
color_cutoff = 'darkorange'
color_new    = 'mediumseagreen'

# ════════════════════════════════════════
# Figure 1: lambda_violin.png (as before)
# ════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# ── Panel 1: overall ──
parts = ax1.violinplot(
    [lambdas_old_all, lambdas_old_hits, lambdas_cutoff],
    positions=[1, 2, 3], showmedians=True, showextrema=True
)
colors1 = ['steelblue', 'seagreen', 'darkorange']
for pc, color in zip(parts['bodies'], colors1):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
parts['cmedians'].set_color('black')
parts['cmedians'].set_linewidth(2)
ax1.axhline(y=1.0, color='gray', linestyle='--', linewidth=1)
ax1.set_xticks([1, 2, 3])
ax1.set_xticklabels(['Original\n(all pairs)', 'Original\n(hit pairs)', 'KING 0.177\n(hit pairs)'])
ax1.set_ylabel('Genomic inflation factor (λGC)')
ax1.set_title('Overall λGC distribution')

# ── Panel 2: per-ancestry (hit pairs) ──
positions_old    = [1, 4, 7]
positions_cutoff = [2, 5, 8]
for i, anc in enumerate(ANCESTRIES):
    if lambdas_by_anc_hits[anc]['old']:
        violin(ax2, lambdas_by_anc_hits[anc]['old'],    positions_old[i],    color_old)
    if lambdas_by_anc_hits[anc]['cutoff']:
        violin(ax2, lambdas_by_anc_hits[anc]['cutoff'], positions_cutoff[i], color_cutoff)
ax2.axhline(y=1.0, color='gray', linestyle='--', linewidth=1)
ax2.set_xticks([1.5, 4.5, 7.5])
ax2.set_xticklabels(ANCESTRIES)
ax2.set_ylabel('Genomic inflation factor (λGC)')
ax2.set_title('λGC by ancestry (hit pairs)')
ax2.legend(handles=[
    mpatches.Patch(facecolor=color_old,    alpha=0.7, label='Original'),
    mpatches.Patch(facecolor=color_cutoff, alpha=0.7, label='KING 0.177'),
])

plt.tight_layout()
plt.savefig('lambda_violin.png', dpi=300)
plt.close()
print("Saved to lambda_violin.png")

# ════════════════════════════════════════
# Figure 2: inflation_comparison.png
# ════════════════════════════════════════
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(18, 6))

# ── Panel 3: overall all pairs — old vs new ──
parts3 = ax3.violinplot(
    [lambdas_old_all, lambdas_new_all],
    positions=[1, 2], showmedians=True, showextrema=True
)
for pc, color in zip(parts3['bodies'], [color_old, color_new]):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
parts3['cmedians'].set_color('black')
parts3['cmedians'].set_linewidth(2)
ax3.axhline(y=1.0, color='gray', linestyle='--', linewidth=1)
ax3.set_xticks([1, 2])
ax3.set_xticklabels(['Original\n(all pairs)', 'New\n(all pairs)'])
ax3.set_ylabel('Genomic inflation factor (λGC)')
ax3.set_title('Overall λGC: Original vs New analysis')
ax3.legend(handles=[
    mpatches.Patch(facecolor=color_old, alpha=0.7, label='Original'),
    mpatches.Patch(facecolor=color_new, alpha=0.7, label='New analysis'),
])

# ── Panel 4: per-ancestry all pairs — old vs new (all 8 ancestries) ──
pos_old4 = [1 + i * 3 for i in range(len(ALL_ANCESTRIES))]
pos_new4 = [2 + i * 3 for i in range(len(ALL_ANCESTRIES))]
for i, anc in enumerate(ALL_ANCESTRIES):
    if lambdas_by_anc_all[anc]['old']:
        violin(ax4, lambdas_by_anc_all[anc]['old'], pos_old4[i], color_old)
    if lambdas_by_anc_all[anc]['new']:
        violin(ax4, lambdas_by_anc_all[anc]['new'], pos_new4[i], color_new)
ax4.axhline(y=1.0, color='gray', linestyle='--', linewidth=1)
ax4.set_xticks([1.5 + i * 3 for i in range(len(ALL_ANCESTRIES))])
ax4.set_xticklabels(ALL_ANCESTRIES)
ax4.set_ylabel('Genomic inflation factor (λGC)')
ax4.set_title('λGC by ancestry (all pairs): Original vs New analysis')
ax4.legend(handles=[
    mpatches.Patch(facecolor=color_old, alpha=0.7, label='Original'),
    mpatches.Patch(facecolor=color_new, alpha=0.7, label='New analysis'),
])

plt.tight_layout()
plt.savefig('inflation_comparison.png', dpi=300)
plt.close()
print("Saved to inflation_comparison.png")
