import pandas as pd
import os

hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS']

hits_list = os.listdir(hit_folder_name)

for hit in [hits_list[0]]:
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/z_scores/results/{hit}'
    os.makedirs(output_folder, exist_ok=True)
    pheno = hit.split('_')[1]
    imp_ancestry = hit.split('_')[0]
    for ancestry in [ancestry_list[0]]:
        file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'

        df = pd.read_csv(file, sep='\t')

        ids = df['ID'].unique().tolist()
        significant_ids = []
        for id in ids:
            add_p_value = df[(df['ID'] == id) & (df['TEST'] == 'ADD')]['P']
            # Filtra per imp_ancestry
            imp_p_value = df[(df['ID'] == id) & (df['TEST'] == imp_ancestry)]['P']

            # Controlla se entrambi i test hanno P < 0.05
            if not add_p_value.empty and not imp_p_value.empty:
                if add_p_value.values[0] < 0.05 and imp_p_value.values[0] < 0.05:
                    significant_ids.append(id)

        # Salva gli SNPs significativi in un file
        output_file = os.path.join(output_folder, f"{hit}_significant_snps_{ancestry}.txt")
        with open(output_file, 'w') as f:
            for snp in significant_ids:
                f.write(f"{snp}\n")

        print(f"Processato {file}. SNPs significativi salvati in {output_file}.")