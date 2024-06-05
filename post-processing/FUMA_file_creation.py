import pandas as pd
import os

type_list = ['no_BMI', 'BMI']

output_folder_generic = '/private/groups/ioannidislab/smeriglio/FUMA/*/'

generic_snps_filename = '/private/groups/ioannidislab/smeriglio/output/*/significant_SNPs_with_single_P_values.txt'

for type in type_list:

    output_folder = output_folder_generic.replace('*', type)
    
    snps_filename = generic_snps_filename.replace('*', type)

    snps = pd.read_csv(snps_filename, sep='\t')

    ancestries = snps['ancestry'].unique()

    phenotypes = snps['phenotype'].unique()

    for ancestry in ancestries:

        for phenotype in phenotypes:

            snps_subset = snps[(snps['ancestry'] == ancestry) & (snps['phenotype'] == phenotype)]

            fuma_snps = snps_subset.drop(columns=['phenotype', 'ancestry', 'start', 'end'])
            fuma_snps.rename(columns={'chr' : 'CHR'}, inplace=True)
            fuma_wind = snps_subset.drop(columns=['pos', 'rfid', 'P', 'phenotype', 'ancestry'])

            output_file_snps = os.path.join(output_folder, 'snps', f'{ancestry}_{phenotype}_snps.txt')
            output_file_wind = os.path.join(output_folder, 'wind', f'{ancestry}_{phenotype}_wind.txt')
            if not snps_subset.empty:
                fuma_snps.to_csv(output_file_snps, sep='\t', index=False)
                fuma_wind.to_csv(output_file_wind, sep='\t', index=False)
