import os
import numpy as np
import pandas as pd

candidates_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
output_file     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_causal_candidates.tsv'
admixture_output_dir = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb'

candidates = pd.read_csv(candidates_file, sep='\t')


# Original admixture mapping OR (LAI window-level ADD test) for each hit's window
_or_cache = {}
def admixture_window_or(hit, window_start):
    key = (hit, window_start)
    if key in _or_cache:
        return _or_cache[key]
    ancestry, pheno, _ = hit.split('_')
    glm_file = os.path.join(admixture_output_dir, f'output_ancestry_{ancestry}', pheno,
                             f'output.{pheno}.glm.logistic.hybrid')
    glm = pd.read_csv(glm_file, sep='\t')
    row = glm[(glm['TEST'] == 'ADD') & (glm['POS'] == window_start)]
    val = row['OR'].iloc[0] if not row.empty else np.nan
    _or_cache[key] = val
    return val


# Concordance with the admixture analysis: the SNP OR direction must match the
# window-level admixture OR direction.
candidates['admixture_OR'] = candidates.apply(
    lambda r: admixture_window_or(r['hit'], r['window_start']), axis=1
)
candidates['concordant_direction'] = (
    np.sign(np.log(candidates['admixture_OR'])) == np.sign(np.log(candidates['OR']))
)

# Keep only candidates concordant with the admixture analysis
concordant = candidates[candidates['concordant_direction']].copy()
print(f'Concordant candidates: {len(concordant)} / {len(candidates)}')

# For each hit, keep the concordant candidate with the largest effect size, measured
# from the OR as |log(OR)| (OR-based distance from null OR=1; identical to |beta|).
concordant['or_effect'] = np.abs(np.log(concordant['OR']))
idx = concordant.groupby('hit')['or_effect'].apply(lambda x: x.idxmax())
top_candidates = concordant.loc[idx, ['hit', 'ID', 'CHROM', 'POS', 'OR', 'beta', 'P_BY', 'LAI_P_BY',
                                       'window_start', 'window_end', 'admixture_OR',
                                       'concordant_direction']]
top_candidates['in_window'] = top_candidates['POS'].between(
    top_candidates['window_start'], top_candidates['window_end']
)
top_candidates = top_candidates.rename(columns={'ID': 'candidate_snp'}).sort_values('hit')

top_candidates.to_csv(output_file, sep='\t', index=False)
print(f'Saved -> {output_file}\n')

print(top_candidates.to_string(index=False))
