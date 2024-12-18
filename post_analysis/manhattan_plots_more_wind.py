import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

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
    output_folder_general = os.path.join(base_folder, 'fine_mapping_ancestries')

    list_folders = os.listdir(output_folder_general)

    df_list = []

    for fold_name in list_folders:

        # Read the file and compute the Bonferroni threshold
        current_out_folder = os.path.join(output_folder_general, fold_name)
        current_file = os.path.join(current_out_folder, f'{fold_name}_output.{fold_name.split("_")[1]}.glm.logistic.hybrid')

        df = pd.read_csv(current_file, sep='\t')
        df_list.append(df)

    index = 0
    for fold_name in list_folders:

        df = df_list[index]
        
        # Find the SNP with the minimum p-value across the entire dataframe
        min_p_value_snp = df.loc[df['P'].idxmin()]

        # Create Manhattan plot for each file
        output_folder = os.path.join(base_folder, 'manhattan_plots_fine_mapping')

        # Create the plot
        plt.figure(figsize=(12, 6))

        # There is only one chromosome per file; extract its label
        chrom = df['#CHROM'].iloc[0]

        plt.scatter(df['POS'], -np.log10(df['P']), 
                    color=chromosome_colors[int(chrom+1) % len(chromosome_colors)], 
                    label=f'Chromosome {chrom}', s=10)

        # Set the x-axis limits based on the data
        plt.xlim(df['POS'].min(), df['POS'].max())

        # Set the x-axis label and ticks
        plt.xlabel(f'Chromosome {chrom}')

        # Plot the Bonferroni significance threshold
        bonferroni_threshold = significance_threshold / len(df)
        plt.axhline(y=-np.log10(bonferroni_threshold), color='r', linestyle='--', label='Bonferroni')

        # Add labels and titles
        plt.title(f'Manhattan Plot {fold_name}')
        plt.ylabel('-log10(p-value)')
        plt.legend(title='Significance Thresholds')

        # Annotate the SNP with the minimum p-value
        plt.annotate(
            f"{min_p_value_snp['ID']}",
            xy=(min_p_value_snp['POS'], -np.log10(min_p_value_snp['P']) + 0.2),
            xytext=(min_p_value_snp['POS'], -np.log10(min_p_value_snp['P']) + 0.2),
            textcoords='data',
            color='black',
            fontsize=9,
            ha='center',
            arrowprops=dict(facecolor='black', arrowstyle="->", connectionstyle="arc3,rad=.5")
        )

        # Save the plot
        plt.tight_layout()
        os.makedirs(output_folder, exist_ok=True)  # Ensure output folder exists
        plt.savefig(os.path.join(output_folder, f'manhattan_plot_{fold_name}.png'))
        plt.close()
        print(f'Manhattan plot of {fold_name} (Chromosome {chrom}) saved to {output_folder}')
        index += 1
