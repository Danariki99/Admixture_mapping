import os
import pandas as pd

# Define the map of significant loci found in AoU
aou_hits = {
    'AFR': {
        'Obesity': {'chrom': 17, 'start': 77484047, 'end': 77663751},
        'obesity': {'chrom': 17, 'start': 77484047, 'end': 77663751},
    },
    'EAS': {
        'diabetes': {'chrom': 6, 'start': 31327241, 'end': 31341458},
        'Diabetes': {'chrom': 6, 'start': 31327241, 'end': 31341458},
    },
    'EUR': {
        'Actinic_keratosis': {'chrom': 5, 'start': 33988445, 'end': 34275649},
        'actinic_keratosis': {'chrom': 5, 'start': 33988445, 'end': 34275649},
        'keratosis': {'chrom': 5, 'start': 33988445, 'end': 34275649},
        'Keratosis': {'chrom': 5, 'start': 33988445, 'end': 34275649},
        'Hypothyroidism': {'chrom': 6, 'start': 31001177, 'end': 31021547},
        'hypothyroidism': {'chrom': 6, 'start': 31001177, 'end': 31021547},
        'diabetes': {'chrom': 6, 'start': 32217361, 'end': 32264645},
        'Diabetes': {'chrom': 6, 'start': 32217361, 'end': 32264645},
        'diabetes': {'chrom': 14, 'start': 101754561, 'end': 101907902},
        'Diabetes': {'chrom': 14, 'start': 101754561, 'end': 101907902},
        'obesity': {'chrom': 11, 'start': 8651131, 'end': 8837745},
        'Obesity': {'chrom': 11, 'start': 8651131, 'end': 8837745},
    },
    'SAS': {
        'spondylosis': {'chrom': 8, 'start': 82999687, 'end': 83271468},
        'Spondylosis': {'chrom': 8, 'start': 82999687, 'end': 83271468},
        'bowel': {'chrom': 3, 'start': 104916249, 'end': 105073444},
        'Bowel': {'chrom': 3, 'start': 104916249, 'end': 105073444},
        'Tobacco': {'chrom': 2, 'start': 204556927, 'end': 204810578},
        'tobacco': {'chrom': 2, 'start': 204556927, 'end': 204810578},
        'diabetes': {'chrom': 16, 'start': 31129895, 'end': 31334236},
        'Diabetes': {'chrom': 16, 'start': 31129895, 'end': 31334236},
    }
}

# Path where UKBB results are stored
result_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb/output_ancestry_*'

# Load the list of phenotypes found in UKBB
df_ukbb = pd.read_excel("../tables_plots/ukbb_v1.xlsx", sheet_name="first_batch")

# Start searching for UKBB replications of AoU findings
ancestries = list(aou_hits.keys())
count = 0

for ancestry in ancestries:
    print(f'Dealing with ancestry {ancestry}')
    result_path_ancestry = result_path.replace('*', ancestry)
    ukbb_phenotypes = df_ukbb['ID2'].tolist()
    ukbb_ids = df_ukbb['ID'].tolist()
    aou_phenotypes = list(aou_hits[ancestry].keys())

    for aou_pheno in aou_phenotypes:
        print(f'Analyzing AoU phenotype {aou_pheno}')
        chrom = aou_hits[ancestry][aou_pheno]['chrom']
        start = aou_hits[ancestry][aou_pheno]['start']
        end = aou_hits[ancestry][aou_pheno]['end']
        compliant_phenotypes = [s for s in ukbb_phenotypes if aou_pheno in s]

        for compliant_pheno in compliant_phenotypes:
            count += 1
            print(f'Analyzing UKBB phenotype {compliant_pheno}')
            index = ukbb_phenotypes.index(compliant_pheno)
            compliant_pheno_id = ukbb_ids[index]
            file_path = os.path.join(result_path_ancestry, compliant_pheno_id, f'output.{compliant_pheno_id}.glm.logistic.hybrid')

            if not os.path.exists(file_path):
                print(f'File not found: {file_path}')
                continue

            df = pd.read_csv(file_path, sep='\t')
            filtered_df = df[(df['#CHROM'] == chrom) & (df['POS'].between(start, end))]
            threshold = 0.05
            significant_rows = filtered_df[filtered_df['P'] <= threshold]

            if not significant_rows.empty:
                print(f'Found a significant row for phenotype {compliant_pheno}')
                print(significant_rows)
            else:
                print(f'No significant row found for phenotype {compliant_pheno}')
            print('\n' * 5)

print(count)