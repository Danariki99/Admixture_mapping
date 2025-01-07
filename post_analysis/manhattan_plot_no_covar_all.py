import os
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, List
import matplotlib.pyplot as plt

def manhattan_plot(
    input_files: List[str],
    colors: List[str],
    significance_threshold: float = 0.05,
    figsize: Optional[Tuple[float, float]] = (20, 10),
    title: Optional[str] = None,
    fontsize: Optional[Dict[str, float]] = None,
    save: Optional[bool] = False,
    output_filename: Optional[str] = "manhattan_plot.png",
):
    if len(input_files) != len(colors):
        raise ValueError("The number of input files must match the number of colors.")

    plt.figure(figsize=figsize)

    # Concatenazione dei file in un unico DataFrame
    combined_df = pd.DataFrame()

    for file, color in zip(input_files, colors):
        df = pd.read_csv(file, sep='\t')
        df['COLOR'] = color
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    # Calcolo della soglia di significatività Bonferroni
    bonferroni_threshold = significance_threshold / len(combined_df)

    # Plot degli SNP
    unique_ids = combined_df['ID'].unique()
    positions = np.arange(len(unique_ids))

    for i, color in enumerate(colors):
        subset = combined_df[combined_df['COLOR'] == color]
        plt.scatter(
            subset['ID'].map(lambda x: np.where(unique_ids == x)[0][0]),
            -np.log10(subset['P']),
            color=color,
            label=f"File {i + 1}"
        )

    # Impostazioni asse X
    plt.xticks(
        ticks=positions[::max(1, len(positions) // 100)],
        labels=unique_ids[::max(1, len(positions) // 100)],
        rotation=90,
        fontsize=fontsize.get('xticks', 8) if fontsize else 8
    )

    # Linea di soglia Bonferroni
    plt.axhline(y=-np.log10(bonferroni_threshold), color='red', linestyle='--', label='Bonferroni threshold')

    # Titolo, etichette e legenda
    if title:
        plt.title(title, fontsize=fontsize.get('title', 20) if fontsize else 20)
    plt.xlabel('SNP IDs', fontsize=fontsize.get('xlabel', 15) if fontsize else 15)
    plt.ylabel('-log10(p-value)', fontsize=fontsize.get('ylabel', 15) if fontsize else 15)
    plt.legend(fontsize=fontsize.get('legend', 12) if fontsize else 12)

    # Salvataggio o visualizzazione del plot
    if save:
        plt.tight_layout()
        plt.savefig(output_filename)
        print(f"Plot saved as: {output_filename}")
    else:
        plt.show()



if __name__ == '__main__':
    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', 
        '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', 
        '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
    ]
    ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS']


    res_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_no_covar_PCA'

    for hit in os.listdir(res_folder):
        current_folder = os.path.join(res_folder, hit)

        input_files = []
        for ancestry in ancestry_list:
            input_files.append(os.path.join(current_folder, f'{hit}_output.{ancestry}.glm.logistic.hybrid'))

        output_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_plots/no_covar/all/manhattan_plot_{hit}.png'

        manhattan_plot(
        input_file=input_files, 
        colors=colors[:len(input_files)],
        significance_threshold=0.05,
        figsize=(45, 8),
        title=f'Manhattan Plot {hit} {hit.split("_")[0]}',
        fontsize={'title': 20, 'xlabel': 8, 'ylabel': 8, 'legend': 8},
        save=True,
        output_filename=output_file
        )

        print(f'Finished {hit}')
        