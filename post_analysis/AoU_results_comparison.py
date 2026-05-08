import os
import pandas as pd

# Loci AoU da cercare in UKBB, con etichette uniche e termini per il match dei fenotipi UKBB
AOU_HITS = {
    'AFR': [
        {'label': 'obesity_chr17', 'terms': ['Obesity', 'obesity'], 'chrom': 17, 'start': 77484047, 'end': 77663751},
    ],
    'EAS': [
        {'label': 'diabetes_chr6', 'terms': ['Diabetes', 'diabetes'], 'chrom': 6, 'start': 31327241, 'end': 31341458},
    ],
    'EUR': [
        {'label': 'actinic_keratosis_chr5', 'terms': ['Actinic_keratosis', 'actinic_keratosis', 'keratosis', 'Keratosis'], 'chrom': 5, 'start': 33988445, 'end': 34275649},
        {'label': 'hypothyroidism_chr6', 'terms': ['Hypothyroidism', 'hypothyroidism'], 'chrom': 6, 'start': 31001177, 'end': 31021547},
        {'label': 'diabetes_chr6', 'terms': ['Diabetes', 'diabetes'], 'chrom': 6, 'start': 32217361, 'end': 32264645},
        {'label': 'diabetes_chr14', 'terms': ['Diabetes', 'diabetes'], 'chrom': 14, 'start': 101754561, 'end': 101907902},
        {'label': 'obesity_chr11', 'terms': ['Obesity', 'obesity'], 'chrom': 11, 'start': 8651131, 'end': 8837745},
    ],
    'SAS': [
        {'label': 'spondylosis_chr8', 'terms': ['Spondylosis', 'spondylosis'], 'chrom': 8, 'start': 82999687, 'end': 83271468},
        {'label': 'bowel_chr3', 'terms': ['Bowel', 'bowel'], 'chrom': 3, 'start': 104916249, 'end': 105073444},
        {'label': 'tobacco_chr2', 'terms': ['Tobacco', 'tobacco'], 'chrom': 2, 'start': 204556927, 'end': 204810578},
        {'label': 'diabetes_chr16', 'terms': ['Diabetes', 'diabetes'], 'chrom': 16, 'start': 31129895, 'end': 31334236},
    ],
}

RESULT_PATH = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb/output_ancestry_*'
P_THRESHOLD = 0.05

df_ukbb = pd.read_excel("../tables_plots/ukbb_v1.xlsx", sheet_name="first_batch")
harmony_dir = os.path.dirname(__file__)
harmonization_table_path = os.path.join(harmony_dir, "aou_ukbb_phenotype_harmonization.tsv")
harmonization_excel_path = os.path.join(harmony_dir, "aou_ukbb_phenotype_harmonization.xlsx")

harmonized_rows = []
total_checks = 0

for ancestry, hits in AOU_HITS.items():
    print(f'Dealing with ancestry {ancestry}')
    result_path_ancestry = RESULT_PATH.replace('*', ancestry)
    if not os.path.isdir(result_path_ancestry):
        print(f'Folder not found: {result_path_ancestry}')
        continue

    ukbb_phenotypes = df_ukbb['ID2'].tolist()
    ukbb_ids = df_ukbb['ID'].tolist()

    for hit in hits:
        chrom = hit['chrom']
        start = hit['start']
        end = hit['end']
        terms = hit['terms']
        aou_label = hit['label']

        compliant_phenotypes = [s for s in ukbb_phenotypes if any(term in s for term in terms)]

        for compliant_pheno in compliant_phenotypes:
            total_checks += 1
            idx = ukbb_phenotypes.index(compliant_pheno)
            compliant_pheno_id = ukbb_ids[idx]

            file_path = os.path.join(result_path_ancestry, compliant_pheno_id, f'output.{compliant_pheno_id}.glm.logistic.hybrid')
            file_exists = os.path.exists(file_path)

            significant_rows = pd.DataFrame()
            hit_info = {
                'ukbb_hit_chrom': None,
                'ukbb_hit_pos': None,
                'ukbb_hit_ref': None,
                'ukbb_hit_alt': None,
                'ukbb_hit_test': None,
                'ukbb_hit_beta': None,
                'ukbb_hit_or': None,
                'ukbb_hit_p': None,
            }

            if not file_exists:
                print(f'File not found: {file_path}')
            else:
                df = pd.read_csv(file_path, sep='\t')
                filtered_df = df[(df['#CHROM'] == chrom) & (df['POS'].between(start, end))]
                significant_rows = filtered_df[filtered_df['P'] <= P_THRESHOLD]

                if not significant_rows.empty:
                    best_row = significant_rows.loc[significant_rows['P'].idxmin()]
                    hit_info = {
                        'ukbb_hit_chrom': best_row.get('#CHROM'),
                        'ukbb_hit_pos': best_row.get('POS'),
                        'ukbb_hit_ref': best_row.get('REF'),
                        'ukbb_hit_alt': best_row.get('ALT'),
                        'ukbb_hit_test': best_row.get('TEST'),
                        'ukbb_hit_beta': best_row.get('BETA'),
                        'ukbb_hit_or': best_row.get('OR'),
                        'ukbb_hit_p': best_row.get('P'),
                    }

            harmonized_rows.append({
                'ancestry': ancestry,
                'aou_pheno_label': aou_label,
                'ukbb_pheno': compliant_pheno,
                'ukbb_id': compliant_pheno_id,
                'aou_chrom': chrom,
                'aou_start': start,
                'aou_end': end,
                'ukbb_file_exists': file_exists,
                'significant_match_found': not significant_rows.empty,
                **hit_info,
            })

print(f'Total checks: {total_checks}')

if harmonized_rows:
    harmonized_df = pd.DataFrame(harmonized_rows)
    harmonized_df.sort_values(by=['ancestry', 'aou_pheno_label', 'ukbb_pheno'], inplace=True, ignore_index=True)
    harmonized_df.to_csv(harmonization_table_path, sep='\t', index=False)
    print(f'Harmonization table saved to: {harmonization_table_path}')

    excel_df = harmonized_df.rename(columns={
        'ancestry': 'Population',
        'aou_pheno_label': 'AoU phenotype label',
        'ukbb_pheno': 'UKBB phenotype',
        'ukbb_id': 'UKBB phenotype ID',
        'aou_chrom': 'AoU hit chrom',
        'aou_start': 'AoU hit start',
        'aou_end': 'AoU hit end',
        'ukbb_hit_chrom': 'UKBB hit chrom',
        'ukbb_hit_pos': 'UKBB hit position',
        'ukbb_hit_ref': 'UKBB REF',
        'ukbb_hit_alt': 'UKBB ALT',
        'ukbb_hit_test': 'TEST (covariate)',
        'ukbb_hit_beta': 'Beta',
        'ukbb_hit_or': 'OR',
        'ukbb_hit_p': 'P-value',
        'significant_match_found': 'Significant match',
        'ukbb_file_exists': 'UKBB file exists',
    })
    excel_df.to_excel(harmonization_excel_path, index=False)
    print(f'Harmonized phenotype pairs saved to: {harmonization_excel_path}')
else:
    print('No AoU ↔ UKBB phenotype matches were recorded; harmonization table not created.')
