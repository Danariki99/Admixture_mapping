import os
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict
import matplotlib.pyplot as plt
import statsmodels.api as sm

def manhattan_plot(
    input_file: str, 
    colors: list,
    significance_threshold: float = 0.05,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    fontsize: Optional[Dict[str, float]] = None,
    save: Optional[bool] = None,
    output_filename: Optional[str] = None,
):
    # Read the input file
    df = pd.read_csv(input_file, sep='\t')

    # Bonferroni threshold
    bonferroni_threshold = significance_threshold / len(df)

    # Create the plot
    plt.figure(figsize=figsize)

    # Display Manhattan plot points for each chromosome
    for i, (chrom, chrom_data) in enumerate(df.groupby('#CHROM')):
        chrom_data['ABS_POS'] = chrom_data['POS']
        plt.scatter(chrom_data['ABS_POS'], -np.log10(chrom_data['P']), 
                    color=colors[int(chrom+1) % len(colors)])

    # X-axis settings
    plt.xticks(df['POS'], df['POS'], rotation=90)

    # Significance thresholds
    plt.axhline(y=-np.log10(bonferroni_threshold), color='r', linestyle='--', label='Bonferroni')

    # Labels and title
    if title:
        plt.title(title, fontsize=fontsize.get('title', 20))
    plt.xlabel(f'Chromosome {chrom}', fontsize=fontsize.get('xlabel', 15))
    plt.ylabel('-log10(p-value)', fontsize=fontsize.get('ylabel', 15))
    plt.legend(fontsize=fontsize.get('legend', 15))

    # Save the plot
    plt.tight_layout()
    if save:
        plt.savefig(output_filename)
    plt.show()

if __name__ == '__main__':
    chromosome_colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', 
        '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', 
        '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
    ]

    res_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA'

    for hit in os.listdir(res_folder):
        current_folder = os.path.join(res_folder, hit)

        input_file = os.path.join(current_folder, f'{hit}_output.{hit.split("_")[0]}.{hit.split("_")[1]}.glm.logistic.hybrid')

        output_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_plots/covar/single/manhattan_plot_{hit}.{hit.split("_")[0]}.png'
    
        manhattan_plot(
        input_file=input_file, 
        colors=chromosome_colors,
        significance_threshold=0.05,  # Default significance threshold
        figsize=(45, 8),  # Optional: Set the figure size
        title=f'Manhattan Plot {hit} {hit.split("_")[0]}',
        fontsize={'title': 20, 'xlabel': 8, 'ylabel': 8, 'legend': 8},
        save=True,
        output_filename=output_file
        )

        print(f'Finished {hit}')
        