import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def manhattan_plot(genetic_file, chr_col, bp_col, p_col, threshold=None, suggestiveline=None, genomewideline=None, title=None):
    # Load data from the genetic file
    data = pd.read_csv(genetic_file, sep='\t', header = 0)

    # update the value of the windows positions
    for chrom in range(2, max(data[chr_col]) + 1):
        previous_chrom_max_pos = data.loc[data[chr_col] == chrom - 1, bp_col].max()
        data.loc[data[chr_col] == chrom, bp_col] += previous_chrom_max_pos

    folder = os.path.dirname(genetic_file)

    filename = os.path.join(folder, 'manhattan_plot.pdf')
    # Calculate -log10(p) for Manhattan plot
    data['-log10p'] = -np.log10(data[p_col])
    chromosome_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', 
    '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', 
    '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
]
    # Create a Manhattan plotù
    pheno = os.path.basename(genetic_file).split('.')[1]
    title = title + ' ' + pheno

    min_pos = data[bp_col].min()
    max_pos = data[bp_col].max()

    plt.figure(figsize=(12, 6))
    # Calculate the chromosome positions and labels
    chrom_positions = []
    chrom_labels = []
    for chrom in data[chr_col].unique():
        if chrom <= 13 or (chrom - 13) % 2 == 0:
            chrom_data = data[data[chr_col] == chrom]
            chrom_positions.append(chrom_data[bp_col].mean())
            chrom_labels.append(chrom)
    # Plot Manhattan points for each chromosome
    for chrom, chrom_data in data.groupby(chr_col):
        plt.scatter(chrom_data[bp_col], -np.log10(chrom_data[p_col]), color=chromosome_colors[chrom-1])
    # Set x-axis ticks to the chromosome positions and labels
    plt.xticks(chrom_positions, chrom_labels)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p)')
    if title:
        plt.title(title)
    if threshold:
        plt.axhline(y=-np.log10(threshold), color='red', linestyle='--', linewidth=1)
    if suggestiveline:
        plt.axhline(y=suggestiveline, color='green', linestyle='--', linewidth=1)
    if genomewideline:
        plt.axhline(y=genomewideline, color='orange', linestyle='--', linewidth=1)
    if filename:
        plt.xlim(min_pos, max_pos)
        plt.savefig(filename)

# Using the function to create a Manhattan plot
manhattan_plot('/private/groups/ioannidislab/smeriglio/output/no_BMI/output_ancestry_EUR/HC1134/output.HC1134.glm.logistic.hybrid', 
               chr_col='#CHROM', bp_col='POS', p_col='P', title='Manhattan Plot')