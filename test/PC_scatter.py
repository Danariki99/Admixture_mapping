import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

filename = '/private/groups/ioannidislab/ukbb24983/sqc/ukb24983_GWAS_covar.phe'

df = pd.read_csv(filename, sep='\t', header=None, low_memory=False)
df.columns = df.iloc[0]
df = df[1:]

# Convert 'PC1' and 'PC2' to numeric, coercing non-numeric values to NaN
df['PC1'] = pd.to_numeric(df['PC1'], errors='coerce')
df['PC2'] = pd.to_numeric(df['PC2'], errors='coerce')

# Similarly for 'Global_PC1' and 'Global_PC2'
df['Global_PC1'] = pd.to_numeric(df['Global_PC1'], errors='coerce')
df['Global_PC2'] = pd.to_numeric(df['Global_PC2'], errors='coerce')

# Drop rows with NaN values in these columns
df = df.dropna(subset=['PC1', 'PC2', 'Global_PC1', 'Global_PC2'])

print(df['population'].value_counts())

# Extract unique categories
categories = df['population'].unique()

# Create a dictionary of dataframes for each population
dfs = {category: df[df['population'] == category] for category in categories}

for i in categories:
    print(f'category: {i} = {len(dfs[i])}')

# Define a list of colors
colors = ['red', 'green', 'blue', 'purple', 'orange']

plt.figure(figsize=(10, 5))
# Iterate over each category
for i, category in enumerate(categories):
    # Select the dataframe for this category
    df_category = dfs[category]
    
    # Create a scatter plot for this category
    plt.scatter(df_category['Global_PC1'], df_category['Global_PC2'], color=colors[i % len(colors)], label=category, alpha=0.5)

plt.title('Scatter plot: Global_PC1 vs Global_PC2')
plt.xlabel('Global_PC1')
plt.ylabel('Global_PC2')
plt.legend(bbox_to_anchor=(0.5, -0.15), loc='upper center', ncol=len(categories))
plt.savefig('scatter_global_correct.png', bbox_inches='tight')

plt.figure(figsize=(10, 5))
# Iterate over each category
for i, category in enumerate(categories):
    # Select the dataframe for this category
    df_category = dfs[category]
    
    # Create a scatter plot for this category
    plt.scatter(df_category['PC1'], df_category['PC2'], color=colors[i % len(colors)], label=category, alpha=0.5)

plt.title('Scatter plot: PC1 vs PC2')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend(bbox_to_anchor=(0.5, -0.15), loc='upper center', ncol=len(categories))
plt.savefig('scatter_PC_correct.png', bbox_inches='tight')