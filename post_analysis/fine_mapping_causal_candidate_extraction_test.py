import os
import sys
import numpy as np
import pandas as pd

# TEST version of fine_mapping_causal_candidate_extraction.py.
# From the concordant fine-mapping candidates, keep for each hit the single
# candidate with the largest effect size (|log(OR)|). Parameterized on result_folder.

if len(sys.argv) < 2:
    print("Usage: python fine_mapping_causal_candidate_extraction_test.py <result_folder>")
    sys.exit(1)

result_folder = sys.argv[1]
cond_dir        = os.path.join(result_folder, 'fine_mapping_conditional_results')
candidates_file = os.path.join(cond_dir, 'fine_mapping_all_candidates.tsv')
output_file     = os.path.join(cond_dir, 'fine_mapping_causal_candidates.tsv')

if not os.path.exists(candidates_file):
    print(f'No candidates file: {candidates_file}. Run fine_mapping_post_processing_test.py first.')
    sys.exit(0)

candidates = pd.read_csv(candidates_file, sep='\t')

# candidates are already concordant (post-processing filtered them); keep concordant
# defensively in case the column is present, then pick the top |log(OR)| per hit.
if 'concordant_direction' in candidates.columns:
    candidates = candidates[candidates['concordant_direction'] == True].copy()   # noqa: E712

candidates['or_effect'] = np.abs(np.log(candidates['OR']))
idx = candidates.groupby('hit')['or_effect'].apply(lambda x: x.idxmax())
keep_cols = [c for c in ['hit', 'ID', 'CHROM', 'POS', 'OR', 'beta', 'P_BY', 'LAI_P_BY',
                         'window_start', 'window_end', 'admixture_OR', 'concordant_direction']
             if c in candidates.columns]
top = candidates.loc[idx, keep_cols].rename(columns={'ID': 'candidate_snp'})
if 'window_start' in top.columns and 'window_end' in top.columns:
    top['in_window'] = top['POS'].between(top['window_start'], top['window_end'])
top = top.sort_values('hit')

top.to_csv(output_file, sep='\t', index=False)
print(f'Saved {len(top)} causal candidates -> {output_file}')
print(top.to_string(index=False))
