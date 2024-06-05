import pandas as pd
import numpy as np
import os
import json

ancestry_list = ["AFR", "AHG", "EAS", "EUR", "NAT", "OCE", "SAS", "WAS"]
#ancestry_list = ["AFR"]
significance_threshold = 0.05

type_list = ['no_BMI', 'BMI']

phe_folder = '/private/groups/ioannidislab/smeriglio/phe_files'
genes_file = '/private/groups/ioannidislab/smeriglio/genes/gene_dict.json'
gene_dict = json.load(open(genes_file))

phe_list = os.listdir(phe_folder)

thresholds_differences_no_BMI = {}
thresholds_differences_BMI = {}

for type in type_list:
    print(type )
    for ancestry in ancestry_list:
        print(ancestry)
        # Define the paths
        data_file = f'/private/groups/ioannidislab/smeriglio/output/{type}/output_ancestry_{ancestry}/P_info_{ancestry}.tsv'

        # Load the data
        data = pd.read_table(data_file, sep="\t")

        # Adjust the significance threshold
        bonferroni_threshold_unadjusted = significance_threshold / len(data)
        bonferroni_threshold_adjusted = significance_threshold / (len(data) * len(phe_list))
        # Create a list to store the significant data for each phenotype
        
        involved_genes_adjusted = []
        involved_genes_unadjusted = []
        flag_dict = {}
        for pheno_file in phe_list:
            pheno = pheno_file.replace('.phe', '')
                    
            pheno_data = data[['#CHROM', 'POS', 'end_POS', pheno]]

            # Filter and sort the data for the current phenotype using the adjusted threshold
            filtered_data_adjusted = pheno_data.loc[pheno_data[pheno] < bonferroni_threshold_adjusted]
            sorted_data_adjusted = filtered_data_adjusted.sort_values(by=pheno)

            # If the DataFrame is not empty, append the involved genes to the adjusted list
            if not sorted_data_adjusted.empty:
                max_row_adjusted = sorted_data_adjusted.loc[sorted_data_adjusted[pheno].idxmin()]
                key_adjusted = f'{max_row_adjusted["#CHROM"].astype(int)}_{max_row_adjusted["POS"].astype(int)}_{max_row_adjusted["end_POS"].astype(int)}'
                #print(key_adjusted)
                for gene in gene_dict[key_adjusted]:
                    involved_genes_adjusted.append(gene['gene_name'])

            # Filter and sort the data for the current phenotype using the unadjusted threshold
            filtered_data_unadjusted = pheno_data.loc[pheno_data[pheno] < bonferroni_threshold_unadjusted]
            sorted_data_unadjusted = filtered_data_unadjusted.sort_values(by=pheno)

            # If the DataFrame is not empty, append the involved genes to the unadjusted list
            if not sorted_data_unadjusted.empty:
                max_row_unadjusted = sorted_data_unadjusted.loc[sorted_data_unadjusted[pheno].idxmin()]
                key_unadjusted = f'{max_row_unadjusted["#CHROM"].astype(int)}_{max_row_unadjusted["POS"].astype(int)}_{max_row_unadjusted["end_POS"].astype(int)}'
                #print(key_unadjusted)
                for gene in gene_dict[key_unadjusted]:
                    involved_genes_unadjusted.append(gene['gene_name'])

        flag_dict['adjusted'] = list(set(involved_genes_adjusted))
        flag_dict['unadjusted'] = list(set(involved_genes_unadjusted))
        flag_dict['lost_genes'] = list(set(involved_genes_unadjusted)- (set(involved_genes_adjusted)))
        #print('Involved genes (adjusted):', flag_dict['adjusted'])
        #print('Involved genes (unadjusted):', flag_dict['unadjusted'])
        #print('lost_genes:', flag_dict['lost_genes'])
        if type == 'no_BMI':
            thresholds_differences_no_BMI[ancestry] = flag_dict
        else:
            thresholds_differences_BMI[ancestry] = flag_dict

# Save the thresholds_differences_no_BMI dictionary to a file
with open('thresholds_differences_no_BMI.json', 'w') as f:
    json.dump(thresholds_differences_no_BMI, f, indent=4)

# Save the thresholds_differences_BMI dictionary to a file
with open('thresholds_differences_BMI.json', 'w') as f:
    json.dump(thresholds_differences_BMI, f, indent=4)


            
        


            