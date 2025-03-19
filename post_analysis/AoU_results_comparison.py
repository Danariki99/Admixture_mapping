import os
import pandas as pd

# Define the map of significant loci found in AoU
aou_hits = {
    'AFR': {
        'Obesity': {'chrom': 17, 'start': 79487964, 'end': 79689827},
        'obesity': {'chrom': 17, 'start': 79487964, 'end': 79689827},
    },
    'EAS': {
        'diabetes': {'chrom': 6, 'start': 31359464, 'end': 31373681},
        'Diabetes': {'chrom': 6, 'start': 31359464, 'end': 31373681},
    },
    'EUR': {
        'Actinic_keratosis': {'chrom': 5, 'start': 33988340, 'end': 34275544},
        'actinic_keratosis': {'chrom': 5, 'start': 33988340, 'end': 34275544},
        'keratosis': {'chrom': 5, 'start': 33988340, 'end': 34275544},
        'Keratosis': {'chrom': 5, 'start': 33988340, 'end': 34275544},
        'Hypothyroidism': {'chrom': 6, 'start': 31033400, 'end': 31053770},
        'hypothyroidism': {'chrom': 6, 'start': 31033400, 'end': 31053770},
        'diabetes': {'chrom': 6, 'start': 32249584, 'end': 32296868},
        'Diabetes': {'chrom': 6, 'start': 32249584, 'end': 32296868},
        'diabetes': {'chrom': 14, 'start': 101288224, 'end': 101441565},
        'Diabetes': {'chrom': 14, 'start': 101288224, 'end': 101441565},
        'obesity': {'chrom': 11, 'start': 8629584, 'end': 8816198},
        'Obesity': {'chrom': 11, 'start': 8629584, 'end': 8816198},
    },
    'SAS': {
        'spondylosis': {'chrom': 8, 'start': 82087452, 'end': 82359233},
        'Spondylosis': {'chrom': 8, 'start': 82087452, 'end': 82359233},
        'bowel': {'chrom': 3, 'start': 105197405, 'end': 105354600},
        'Bowel': {'chrom': 3, 'start': 105197405, 'end': 105354600},
        'Tobacco': {'chrom': 2, 'start': 203692204, 'end': 203945855},
        'tobacco': {'chrom': 2, 'start': 203692204, 'end': 203945855},
        'diabetes': {'chrom': 16, 'start': 31118574, 'end': 31322915},
        'Diabetes': {'chrom': 16, 'start': 31118574, 'end': 31322915},
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