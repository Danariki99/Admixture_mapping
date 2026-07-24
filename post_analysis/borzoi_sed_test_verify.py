"""
Verify the Borzoi SED TEST run: for the 5 test SNPs, check that sed_test/f3c0 now
has NON-ZERO SED and reproduces the old (valid) sed_all/f3c0 values.
Run in the borzoi_env AFTER borzoi_sed_test.sbatch finishes.
"""
import h5py, numpy as np

B = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results'
TEST = f'{B}/sed_test/f3c0/sed.h5'
OLD  = f'{B}/sed_all/f3c0/sed.h5'
SNPS = ['rs61737338', 'rs3177928', 'rs28732226', 'rs6906021', 'rs9274569']

def per_snp_gene_max(path):
    with h5py.File(path, 'r') as f:
        snp = [x.decode() for x in f['snp'][()]]; si = f['si'][()]
        gene = [x.decode() for x in f['gene'][()]]
        out = {}
        for p in range(len(si)):
            s = snp[si[p]]
            if s in SNPS:
                out.setdefault(s, {})[gene[p]] = f['SED'][p]
    return out

test = per_snp_gene_max(TEST)
old  = per_snp_gene_max(OLD)

print(f"{'SNP':<12}{'genes':<8}{'test max|SED|':<16}{'old max|SED|':<15}{'corr(top gene)'}")
for s in SNPS:
    if s not in test:
        print(f'{s:<12} ASSENTE nel test'); continue
    tmax = max(np.max(np.abs(v)) for v in test[s].values())
    omax = max(np.max(np.abs(v)) for v in old.get(s, {}).values()) if s in old else float('nan')
    # correlation on the strongest-effect gene (by old)
    corr = float('nan')
    if s in old:
        topg = max(old[s], key=lambda g: np.max(np.abs(old[s][g])))
        if topg in test[s]:
            a = old[s][topg].astype(float); w = test[s][topg].astype(float)
            if a.std() > 0 and w.std() > 0:
                corr = np.corrcoef(a, w)[0, 1]
    flag = ' <-- ZERO!' if tmax < 1e-6 else ''
    print(f'{s:<12}{len(test[s]):<8}{tmax:<16.2f}{omax:<15.2f}{corr:.4f}{flag}')

print('\nOK se: test max|SED| NON e\' zero e corr ~1 con il vecchio.')
