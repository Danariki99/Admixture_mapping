import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ancestry_list = ["AFR", "AHG", "EAS", "EUR", "NAT", "OCE", "SAS", "WAS"]

#ancestry_list = ['OCE']
for ancestry in ancestry_list:

    print(ancestry)
    
    # Define the paths
    phe_folder = '/private/groups/ioannidislab/smeriglio/phe_files'
    general_file = f'/private/groups/ioannidislab/smeriglio/output/no_BMI/output_ancestry_{ancestry}/*/output.*.glm.logistic.hybrid'
    window_pos_file = '/private/groups/ioannidislab/smeriglio/windows_pos/positions.csv'

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

    # Filter and sort the data for the first phenotype
    first_filtered_data = data.loc[data[first_pheno] < 0.000005]
    sorted_data = first_filtered_data.sort_values(by=first_pheno)

    # Store the significant data
    significant_list = [sorted_data]

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
        filtered_data = df.loc[df[pheno] < 0.000005]
        sorted_data = filtered_data.sort_values(by=pheno)

        # Store the significant data
        significant_list.append(sorted_data)

        # Merge the data for the current phenotype into the main dataframe
        data = pd.merge(data, df, on=['#CHROM', 'POS', 'ABS_POS'])
    print(data)

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
            plt.annotate(f'{pheno}_{max_row["POS"].astype(int)}', (max_row['ABS_POS'], -np.log10(max_row[pheno])))

    # Add a horizontal line at -log10(0.000005)
    plt.axhline(y=-np.log10(0.000005), color='r', linestyle='--')

    # Set x-axis ticks to the chromosome positions and labels
    plt.xticks(chrom_positions, chrom_labels)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p)')
    plt.title('Manhattan Plot_' + ancestry)
    plt.savefig('../plots/manhattan_plot_' + ancestry)
