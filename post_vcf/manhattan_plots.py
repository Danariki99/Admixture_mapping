import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    # Check if the dataset argument is provided
    if len(sys.argv) != 2:
        print("Usage: python post_processing.py <dataset>")
        sys.exit(1)

    # The dataset variable is taken from the command line argument
    dataset = sys.argv[1]

    # Define the colors for each chromosome
    chromosome_colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', 
        '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', 
        '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
    ]

    significance_threshold = 0.05
    
    base_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/{dataset}'
    output_folder_general = os.path.join(base_folder, 'output')

    list_folders = os.listdir(output_folder_general)

    df_list = []

    max_distance = 0  

    for fold_name in list_folders:

        # Read the file and compute the Bonferroni threshold
        current_out_folder = os.path.join(output_folder_general, fold_name)
        current_file = os.path.join(current_out_folder, f'{fold_name}_output.{fold_name.split("_")[1]}.glm.logistic.hybrid')

        df = pd.read_csv(current_file, sep='\t')

        # Group by chromosome and calculate the distance (max_pos - min_pos) for each chromosome
        for chrom, chrom_data in df.groupby('#CHROM'):
            chrom_max_pos = chrom_data['POS'].max()  # Maximum position in this chromosome
            chrom_min_pos = chrom_data['POS'].min()  # Minimum position in this chromosome
            distance = chrom_max_pos - chrom_min_pos  # Distance between max and min positions

            # Update the global maximum distance if the current distance is greater
            if distance > max_distance:
                max_distance = distance
        df_list.append(df)

    # Calculate the absolute position for each SNP
    for df in df_list:
        df['ABS_POS'] = df['POS'] + max_distance * df['#CHROM']

    index = 0
    for fold_name in list_folders:

        df = df_list[index]
        # Bonferroni correction
        bonferroni_threshold = significance_threshold / len(df)

        # Create Manhattan plot for each file
        output_folder = os.path.join(base_folder, 'manhattan_plots')

        # Create the plot
        plt.figure(figsize=(12, 6))

        # Plot Manhattan points for each chromosome, assigning a unique color to each chromosome
        for i, (chrom, chrom_data) in enumerate(df.groupby('#CHROM')):
            plt.scatter(chrom_data['ABS_POS'], -np.log10(chrom_data['P']), 
                        color=chromosome_colors[i % len(chromosome_colors)], label=f'Chromosome {chrom}', s=10)

        # Plot the Bonferroni significance threshold
        plt.axhline(y=-np.log10(bonferroni_threshold), color='r', linestyle='--')

        # Add labels and titles
        plt.title(f'Manhattan Plot {fold_name}')
        plt.xlabel('Chromosome Position (ABS_POS)')
        plt.ylabel('-log10(p-value)')
        plt.legend(title='Chromosomes')

        # Save the plot
        plt.tight_layout()
        os.makedirs(output_folder, exist_ok=True)  # Ensure output folder exists
        plt.savefig(os.path.join(output_folder, f'manhattan_plot_{fold_name}.png'))
        plt.close()
        print(f'Manhattan plot of {fold_name} saved to {output_folder}')
        index += 1











