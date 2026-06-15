import os
import numpy as np
import pandas as pd

candidates_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
output_file     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_causal_candidates.tsv'
admixture_output_dir = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb'

candidates = pd.read_csv(candidates_file, sep='\t')

# For each hit, keep the candidate with the largest effect size, measured from the
# OR as |log(OR)| (the OR-based distance from the null OR=1; identical to |beta|).
candidates['or_effect'] = np.abs(np.log(candidates['OR']))
idx = candidates.groupby('hit')['or_effect'].apply(lambda x: x.idxmax())
top_candidates = candidates.loc[idx, ['hit', 'ID', 'CHROM', 'POS', 'OR', 'beta', 'P_BY', 'LAI_P_BY',
                                       'window_start', 'window_end']]
top_candidates['in_window'] = top_candidates['POS'].between(
    top_candidates['window_start'], top_candidates['window_end']
)
top_candidates = top_candidates.rename(columns={'ID': 'candidate_snp'}).sort_values('hit')


# Original admixture mapping OR (LAI window-level ADD test) for each hit's window
def admixture_window_or(hit, window_start):
    ancestry, pheno, _ = hit.split('_')
    glm_file = os.path.join(admixture_output_dir, f'output_ancestry_{ancestry}', pheno,
                             f'output.{pheno}.glm.logistic.hybrid')
    glm = pd.read_csv(glm_file, sep='\t')
    row = glm[(glm['TEST'] == 'ADD') & (glm['POS'] == window_start)]
    if row.empty:
        return np.nan
    return row['OR'].iloc[0]


top_candidates['admixture_OR'] = top_candidates.apply(
    lambda r: admixture_window_or(r['hit'], r['window_start']), axis=1
)
top_candidates['concordant_direction'] = (
    np.sign(np.log(top_candidates['admixture_OR'])) == np.sign(np.log(top_candidates['OR']))
)

top_candidates.to_csv(output_file, sep='\t', index=False)
print(f'Saved -> {output_file}\n')

print(top_candidates.to_string(index=False))
