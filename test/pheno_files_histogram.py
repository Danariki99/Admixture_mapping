import os
import pandas as pd
import matplotlib.pyplot as plt

pheno_folder = '/private/groups/ioannidislab/smeriglio/phe_files/'
output_folder = '/private/groups/ioannidislab/smeriglio/phenotype_histograms'

keep_files = '/private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt'

keep = pd.read_csv(keep_files)

keep['#IID'] = keep['#IID'].astype(int)

pheno_files_list = os.listdir(pheno_folder)

# Create two figures and sets of subplots
fig_cases, ax_cases = plt.subplots(figsize=(10, 6))
fig_controls, ax_controls = plt.subplots(figsize=(10, 6))

total_rows = len(keep)

counter = 0

for pheno_file in pheno_files_list:
    
    pheno = pheno_file.replace('.phe', '')

    current_pheno_path = os.path.join(pheno_folder, pheno_file)

    try:
        df = pd.read_csv(current_pheno_path, sep=' ')
        df['#IID'] = df['#IID'].astype(int)
    except:
        # Read the DataFrame without the first line (header)
        df = pd.read_csv(current_pheno_path, sep='\t', header=None, skiprows=1)

        # Manually set the column names
        df.columns = ['#IID', pheno]

        df['#IID'] = df['#IID'].astype(int) 

    df = pd.merge(df, keep, on='#IID')

    df['#IID'] = df['#IID'].astype(str).str.strip()
    if len(df) != 441814:
        print(f'for the pheno {pheno} the number of rows is {len(df)}')

    counts = df[pheno].value_counts()

    # Normalize the counts and create a bar for the count of '2' in the 'HC1230' column with a smaller width
    ax_cases.bar(pheno, counts[2]/total_rows, width=0.8)
    ax_controls.bar(pheno, counts[1]/total_rows, width=0.8)

    counter +=1

# Rotate the x-axis labels and reduce the font size
plt.setp(ax_cases.get_xticklabels(), rotation=90, fontsize=6)
plt.setp(ax_controls.get_xticklabels(), rotation=90, fontsize=6)

# Set the title of the plots
ax_cases.set_title('Cases')
ax_controls.set_title('Controls')

# Save the figures
figure_name_cases = os.path.join(output_folder, 'cases.png')
figure_name_controls = os.path.join(output_folder, 'controls.png')
fig_cases.savefig(figure_name_cases)
fig_controls.savefig(figure_name_controls)

plt.close('all')

print(counter)

# these are the phenotype that i wasn't able to read in this way: [FH1001 FH1220 FH1019 FH1065 FH1286 HC1003 FH1113 FH1263 HC1002 FH1002 FH1044 FH1262]







