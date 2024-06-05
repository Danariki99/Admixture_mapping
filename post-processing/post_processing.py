import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ancestry_list = ["AFR", "AHG", "EAS", "EUR", "NAT", "OCE", "SAS", "WAS"]

#OCE has been removed for lack of data
ancestry_list = ["AFR", "AHG", "EAS", "EUR", "NAT", "SAS", "WAS"]
#ancestry_list = ["AFR"]

significance_threshold = 0.05

type_list = ['no_BMI', 'BMI']
#type_list = ['no_BMI']

for type in type_list:
    print(type)

    significant_df = pd.DataFrame(columns=['#CHROM', 'POS', 'end_POS', 'ABS_POS', 'P', 'Phenotype', 'Ancestry'])

    #ancestry_list = ['OCE']
    for ancestry in ancestry_list:
        print(ancestry)
        
        # Define the paths
        phe_folder = '/private/groups/ioannidislab/smeriglio/phe_files'
        general_file = f'/private/groups/ioannidislab/smeriglio/output/{type}/output_ancestry_{ancestry}/*/output.*.glm.logistic.hybrid'
        window_pos_file = '/private/groups/ioannidislab/smeriglio/windows_pos/positions.csv'
        output_file = f'/private/groups/ioannidislab/smeriglio/output/{type}/output_ancestry_{ancestry}/P_info_{ancestry}.tsv'

        # Get the list of phenotype files
        pheno_list = os.listdir(phe_folder)

        # Process the first phenotype
        first_pheno = pheno_list[0].replace('.phe', '')
        first_file = general_file.replace('*', first_pheno)

        # Load the data for the first phenotype
        data = pd.read_table(first_file, sep="\t")
        data = data[['#CHROM', 'POS', 'P']]
        data = data.rename(columns={'P':first_pheno})

        # Calculate the maximum position for each chromosome
        max_pos = data.groupby('#CHROM')['POS'].max().cumsum()

        # Shift the max_pos Series so that the first value is 0
        max_pos = max_pos.shift(fill_value=0)

        # Calculate ABS_POS directly using the max_pos map
        data['ABS_POS'] = data['POS'] + data['#CHROM'].map(max_pos)

        # Read the window positions file
        window_pos = pd.read_csv(window_pos_file, sep='\t')

        # Merge the data with the window positions
        data = pd.merge(data, window_pos, on=['#CHROM', 'POS'], how='left')

        # Reorder the columns
        data = data[['#CHROM', 'POS', 'end_POS', 'ABS_POS', first_pheno]]

        # Calculate the mean absolute position for each chromosome
        chrom_positions = data.groupby('#CHROM')['ABS_POS'].mean().tolist()

        # Get the chromosome labels
        chrom_labels = sorted(data['#CHROM'].unique())

        # compute bonferroni threshold
        #bonferroni_threshold = significance_threshold / (len(data) * len(pheno_list ) * len(ancestry_list))
        bonferroni_threshold = significance_threshold / len(data)

        # Filter and sort the data for the first phenotype
        first_filtered_data = data.loc[data[first_pheno] < bonferroni_threshold]
        sorted_data = first_filtered_data.sort_values(by=first_pheno)

        # Store the significant data
        significant_list = [sorted_data]

    

        #bonferroni_list = [bonferroni_threshold]

        # Process the remaining phenotypes
        for phe_file in pheno_list[1:]:
            pheno = phe_file.replace('.phe', '')
            current_file = general_file.replace('*', pheno)

            # Load the data for the current phenotype
            df = pd.read_table(current_file, sep="\t")
            df = df[['#CHROM', 'POS', 'P']]
            df = df.rename(columns={'P':pheno})

            # Calculate ABS_POS directly using the max_pos map
            df['ABS_POS'] = df['POS'] + df['#CHROM'].map(max_pos)

            # Filter and sort the data for the current phenotype
            filtered_data = df.loc[df[pheno] < bonferroni_threshold]
            filtered_data = pd.merge(filtered_data, window_pos, on=['#CHROM', 'POS'], how='left')
            sorted_data = filtered_data.sort_values(by=pheno)

            # Store the significant data
            significant_list.append(sorted_data)

            # Merge the data for the current phenotype into the main dataframe
            data = pd.merge(data, df, on=['#CHROM', 'POS', 'ABS_POS'])
        #print(data)

        # save the data as a csv file
        data.to_csv(output_file, index=False, sep = '\t')


        #print(bonferroni_list)

        # Define the colors for each chromosome
        chromosome_colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
            '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', 
            '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', 
            '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
        ]

        # Create the plot
        plt.figure(figsize=(12, 6))

        # Iterate over the first two phenotypes
        for i, phe_file in enumerate(pheno_list):
            pheno = phe_file.replace('.phe', '')

            # Plot Manhattan points for each chromosome
            for chrom, chrom_data in data.groupby('#CHROM'):
                plt.scatter(chrom_data['ABS_POS'], -np.log10(chrom_data[pheno]), color=chromosome_colors[chrom-1], label=pheno)

            # Annotate significant points
            significant_data = significant_list[i]
            if not significant_data.empty:
                max_row = significant_data.loc[significant_data[pheno].idxmin()]
                offset_x = np.random.uniform(-1e6, 1e6)  # Adjust these values as needed
                offset_y = np.random.uniform(-0.5, 0.5)  # Adjust these values as needed
                plt.annotate(f'{pheno}_chr{max_row["#CHROM"].astype(int)}_{max_row["POS"].astype(int)}_{max_row["end_POS"].astype(int)}', (max_row['ABS_POS'] + offset_x, -np.log10(max_row[pheno]) + offset_y))        # Add a horizontal line at -log10(significance_threshold)

                # Rename the phenotype column to 'P'
                significant_data = significant_data.rename(columns={pheno: 'P'})

                # Add the phenotype and ancestry columns
                significant_data['Phenotype'] = pheno
                significant_data['Ancestry'] = ancestry

                # Append the data to the main DataFrame
                significant_df = pd.concat([significant_df, significant_data])

        plt.axhline(y=-np.log10(bonferroni_threshold), color='r', linestyle='--')

        # Set x-axis ticks to the chromosome positions and labels
        plt.xticks(chrom_positions, chrom_labels)
        plt.xlabel('Chromosome')
        plt.ylabel('-log10(p)')
        plt.title('Manhattan Plot_' + ancestry)
        plt.savefig(f'../plots/unadjusted_threshold/{type}/manhattan_plot_' + ancestry)

    print(significant_df)
    #significant_df.to_csv(f'/private/groups/ioannidislab/smeriglio/output/{type}/significant_positions.tsv', sep='\t', index=False)
    
