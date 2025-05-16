import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO
from adjustText import adjust_text
import requests


def positions_extraction(input_file, output_folder):
    # Define the command
    command = f"awk -F'\t' '!/^##/ {{print $1\"\t\"$2}}' {input_file}"

    # Run the command
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    output, error = process.communicate()

    data = StringIO(output.decode())
    df = pd.read_csv(data, sep="\t")

    # Split the 'POS' column into 'start_pos' and 'end_pos'
    df[['POS', 'end_POS']] = df['POS'].str.split('_', expand=True)

    # Create output file path
    output_file = os.path.join(output_folder, 'positions.csv')
    df.to_csv(output_file, sep='\t', index=False)

    return output_file

def result_analysis(ancestry_list, phe_folder, general_file_ini, window_pos_file, general_output_file, plot_output_folder, general_output_folder):
    

    significance_threshold = 0.05

    #read excel file for ukbb (modify for allofus)
    excel_df = pd.read_excel(
        '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx', 
        sheet_name='first_batch',
        usecols="B:C"
    )

    significant_df = pd.DataFrame(columns=['#CHROM', 'POS', 'end_POS', 'ABS_POS', 'P', 'Phenotype', 'Ancestry'])

    #ancestry_list = ['OCE']
    for ancestry in ancestry_list:
        print(ancestry)
        
        # Define the paths
        general_file = general_file_ini.replace('#', ancestry)
        output_file = general_output_file.replace('#', ancestry)

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
        bonferroni_threshold = significance_threshold / (len(data) * len(pheno_list) * 7)
        local_bonferroni_threshold = significance_threshold / len(data)
        
        # Filter and sort the data for the first phenotype
        first_filtered_data = data.loc[data[first_pheno] < bonferroni_threshold]
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
        texts = []
        for i, phe_file in enumerate(pheno_list):
            pheno = phe_file.replace('.phe', '')

            # Plot Manhattan points for each chromosome
            for chrom, chrom_data in data.groupby('#CHROM'):
                plt.scatter(chrom_data['ABS_POS'], -np.log10(chrom_data[pheno]), color=chromosome_colors[chrom-1], label=pheno, s = 7)

            # Annotate significant points
            significant_data = significant_list[i]
            if not significant_data.empty:
                max_row = significant_data.loc[significant_data[pheno].idxmin()]
                if ancestry == 'SAS':
                    value = 0.2
                else:
                    value = 0.6

                offset_x = np.random.uniform(-1e9, 1e9)
                offset_y = np.random.uniform(-value, value)
                name = fetch_cytoband(f"chr{int(max_row['#CHROM'])}", int(max_row["POS"]), int(max_row["end_POS"]))

                # manually select the non working genes
                if name == 'Timeout':
                    if f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr6:31346445-31377047':
                        name = '6p21.33'
                    elif f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr8:124070432-124092625':
                        name = '8q24.13'
                    elif f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr6:31905130-32007956':
                        name = '6p21.33'
                    elif f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr6:32207393-32288190':
                        name = '6p21.32'
                    elif f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr10:116036889-116139029':
                        name = '10q25.3'
                    elif f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr6:31428169-31435326':
                        name = '6p21.33'
                    elif f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr9:85752837-85810910':
                        name = '9q21.32'
                    elif f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}' == 'chr17:1820750-1925859':
                        name = '17p13.3'
                if name != 'Timeout':
                    annotation_text = (
                        f'{excel_df.loc[excel_df["ID"] == pheno, "ID2"].iloc[0].replace("_", " ")}\n'
                        f'{name}'
                    )
                else:
                    annotation_text = (
                        f'{excel_df.loc[excel_df["ID"] == pheno, "ID2"].iloc[0].replace("_", " ")}\n'
                        f'chr{int(max_row["#CHROM"])}:{int(max_row["POS"])}-{int(max_row["end_POS"])}'
                    )

                text = plt.annotate(
                    annotation_text,
                    (max_row['ABS_POS'] + offset_x, -np.log10(max_row[pheno]) + offset_y)
                )
                texts.append(text)
                

                # Rename the phenotype column to 'P'
                significant_data = significant_data.rename(columns={pheno: 'P'})

                # Add the phenotype and ancestry columns
                significant_data['Phenotype'] = pheno
                significant_data['Ancestry'] = ancestry

                # Append the data to the main DataFrame
                significant_df = significant_df.dropna()
                significant_df = pd.concat([significant_df, significant_data])

        adjust_text(texts)
        line1 = plt.axhline(y=-np.log10(bonferroni_threshold), color='r', linestyle='--', label='Bonferroni threshold')
        line2 = plt.axhline(y=-np.log10(local_bonferroni_threshold), color='b', linestyle='--', label='Local Bonferroni threshold')

        # Set x-axis ticks to the chromosome positions and labels
        plt.xticks(chrom_positions, chrom_labels)
        plt.xlabel('Chromosome')
        plt.ylabel('-log10(p)')
        if ancestry == 'NAT':
            ancestry = 'AMR'
        plt.title('Manhattan Plot ' + ancestry)
        plt.legend(handles=[line1, line2], loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
        plt.savefig(os.path.join(plot_output_folder, f'manhattan_plot_{ancestry}.png'), bbox_inches='tight')
    significant_file = os.path.join(general_output_folder, 'significant_positions.tsv')
    if not significant_df.empty:
        significant_df.to_csv(significant_file, sep='\t', index=False)
        return significant_file
    else:
        return None


def SNPs_extraction(input_file, output_dir):
    # Call the R script using subprocess
    current_dir = os.getcwd()
    
    result = subprocess.run(['Rscript', os.path.join(current_dir, 'SNPs_gene_extraction.R'), '-i', input_file, '-o', output_dir], capture_output=True, text=True)

    # Split the stdout into lines and get the last line (the output file path)
    output_lines = result.stdout.strip().split('\n')
    output_file_path = output_lines[-1] if output_lines else ""

    return output_file_path


# Create a function to find the closest SNP to the middle of a given window
def find_snps_in_window(window, snps_df):
    snps_in_window = snps_df[snps_df['pos'].between(window['POS'], window['end_POS'])].copy()
    if not snps_in_window.empty:
        snps_in_window['Phenotype'] = window['Phenotype']
        snps_in_window['Ancestry'] = window['Ancestry']
        snps_in_window['window_start'] = window['POS']
        snps_in_window['window_end'] = window['end_POS']
        snps_in_window['P'] = window['P']
        return snps_in_window[['CHR', 'pos', 'rfid', 'P', 'Phenotype', 'Ancestry', 'window_start', 'window_end']]
    return pd.DataFrame(columns=['CHR', 'pos', 'rfid', 'P', 'Phenotype', 'Ancestry', 'window_start', 'window_end'])

def associate_SNPs_to_windows(snps_file, window_file, output_folder):
    # Load the snps file
    snps_df = pd.read_csv(snps_file, sep='\t')
    snps_df.rename(columns={'chr_name': 'CHR'}, inplace=True)
    snps_df.rename(columns={'refsnp_id': 'rfid'}, inplace=True)
    snps_df.rename(columns={'chrom_start': 'pos'}, inplace=True)
    snps_df = snps_df.drop(columns=['allele'])

    window_df = pd.read_csv(window_file, sep='\t')

    # Convert the windows to intervals and calculate the mid point
    window_df['window'] = pd.IntervalIndex.from_arrays(window_df['POS'], window_df['end_POS'], closed='both')
    window_df['mid'] = window_df['window'].map(lambda x: x.mid)

    all_snps = pd.DataFrame()
    for _, window in window_df.iterrows():
        snps_in_window = find_snps_in_window(window, snps_df)
        all_snps = pd.concat([all_snps, snps_in_window])

    all_snps.columns = ['chr', 'pos', 'rfid', 'P', 'phenotype', 'ancestry', 'start', 'end']

    # Save the new file
    output_file = os.path.join(output_folder, "significant_SNPs_with_P_values.txt")
    all_snps.to_csv(output_file, sep='\t', index=False)
    return output_file

def FUMA_files_creation(snps_filename, output_folder):

    # create the necessary folders
    output_folder_snps = os.path.join(output_folder, 'snps')
    if not os.path.exists(output_folder_snps):
        os.makedirs(output_folder_snps)

    output_folder_wind = os.path.join(output_folder, 'wind')
    if not os.path.exists(output_folder_wind):
        os.makedirs(output_folder_wind)

    snps = pd.read_csv(snps_filename, sep='\t')

    ancestries = snps['ancestry'].unique()

    phenotypes = snps['phenotype'].unique()

    for ancestry in ancestries:

        for phenotype in phenotypes:

            snps_subset = snps[(snps['ancestry'] == ancestry) & (snps['phenotype'] == phenotype)]

            fuma_snps = snps_subset.drop(columns=['phenotype', 'ancestry', 'start', 'end'])
            fuma_snps.rename(columns={'chr' : 'CHR'}, inplace=True)
            fuma_wind = snps_subset.drop(columns=['pos', 'rfid', 'P', 'phenotype', 'ancestry'])
            fuma_wind = fuma_wind.drop_duplicates()

            output_file_snps = os.path.join(output_folder_snps, f'{ancestry}_{phenotype}_snps.txt')
            output_file_wind = os.path.join(output_folder_wind, f'{ancestry}_{phenotype}_wind.txt')
            if not snps_subset.empty:
                fuma_snps.to_csv(output_file_snps, sep='\t', index=False)
                fuma_wind.to_csv(output_file_wind, sep='\t', index=False)
            
    return output_folder_snps, output_folder_wind

def fetch_cytoband(chromosome, start, end, genome="hg19"):
    url = "https://genome.ucsc.edu/cgi-bin/hgTables"
    params = {
        "db": genome,
        "hgta_group": "allTracks",
        "hgta_track": "cytoBand",
        "hgta_table": "cytoBand",
        "hgta_regionType": "range",
        "position": f"{chromosome}:{start}-{end}",
        "hgta_outputType": "primaryTable",
        "boolshad.sendToGalaxy": "0",
        "boolshad.sendToGreat": "0",
        "hgta_doTopSubmit": "get output",
    }

    try:
        response = requests.post(url, data=params, timeout=60)
        response.raise_for_status()  # Gestisce errori HTTP
    except requests.exceptions.ReadTimeout:
        print(f"Timeout: il server non ha risposto entro il tempo previsto.")
        return "Timeout"
    except requests.exceptions.RequestException as e:
        print(f"Errore nella richiesta: {e}")
        return "Errore"
    
    # Analizza la risposta
    response_text = response.text.strip()
    if not response_text:
        return "Unknown"
    
    lines = response_text.split("\n")
    fields = lines[1].split("\t")
    return f"{fields[0].replace('chr', '')}{fields[3]}"

    
