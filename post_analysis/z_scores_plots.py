import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ratios_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/z_scores/results/'
ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS']

for hit in os.listdir(ratios_folder_name):
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/z_scores/plots/{hit}'
    os.makedirs(output_folder, exist_ok=True)

    for ancestry in ancestry_list:
        input_file = f'{ratios_folder_name}{hit}/{hit}_output.{ancestry}.z_scores.txt'
        plot_name = f'{output_folder}/{hit}_output.{ancestry}.z_scores.png'

        df = pd.read_csv(input_file, sep='\t')

        plt.figure(figsize=(45, 6))

        plt.scatter(df['ID'], df['Z_STAT'], alpha=0.6, color='blue', s=30)

        plt.axhline(1, color='red', linestyle='--', label='Z_STAT = 1')
        plt.xlabel('ID')
        plt.ylabel('Z_STAT')
        plt.title(f'Scatter Plot for {hit} - {ancestry}')
        plt.legend()

        x_ticks = np.arange(0, len(df), 1)
        plt.xticks(x_ticks, df['ID'], rotation=90, fontsize=8)

        plt.xlim(0, x_ticks[-1] + 1)

        plt.tight_layout()
        plt.savefig(plot_name)
        plt.close()

        print(f"Plot saved as: {plot_name}")
