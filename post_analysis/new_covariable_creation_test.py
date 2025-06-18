import sys
import subprocess
import os
import pandas as pd
import polars as pl
import snputils as su
import numpy as np


if __name__ == '__main__':
    # Check if the dataset argument is provided
    if len(sys.argv) != 3:
        print("Usage: python post_processing.py <wind_filename>")
        sys.exit(1)
    
    result_folder = sys.argv[1]
    data_folder = sys.argv[2] 


    # Use the dataset variable to construct file paths
    output_folder = os.path.join(result_folder, 'wind_covar_files_processed')
    os.makedirs(output_folder, exist_ok=True)
    wind_folder = os.path.join(result_folder, 'FUMA/wind')
    tmp_folder = os.path.join(result_folder, 'tmp')
    generic_msp_file = os.path.join(result_folder, 'msp_files/chr*.msp')
    old_covar_file = os.path.join(data_folder, 'input.covar')

    # Load old covariate file
    old_covar_df = pd.read_csv(old_covar_file, sep='\t')

    # The dataset variable is taken from the command line argument
    wind_files = os.listdir(wind_folder)
    for wind_filename in wind_files:

        ancestry = wind_filename.split('_')[0]
        pheno = wind_filename.split('_')[1]
        input_file = os.path.join(result_folder, f'vcf_files/ancestry_{ancestry}.vcf')
        wind_file = os.path.join(wind_folder, wind_filename)
        wind_df_1 = pd.read_csv(wind_file, sep='\t')
        chroms = wind_df_1['chr'].unique()

        for chr in chroms:
            msp_file = generic_msp_file.replace('*', str(chr))
            output_file_1 = os.path.join(tmp_folder, f'{pheno}_{ancestry}_{chr}_1.tsv')

            # Awk command to skip lines and extract CHROM and POS columns
            print(f"Running awk command for chromosome {chr}")
            awk_command = f"awk '!/^##/ && $1 == \"{chr}\" {{print $1, $2}}' {input_file} > {output_file_1}"
            subprocess.run(awk_command, shell=True, check=True)

            # Process the output_file_1 to extract ranges
            df = pd.read_csv(output_file_1, sep=' ', header=None, names=['chr', 'pos'])
            df[['start', 'end']] = df['pos'].str.split('_', expand=True)
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            df = df.drop(columns=['pos'])

            wind_df = wind_df_1.loc[wind_df_1['chr'] == chr]
            n_wind = 3

            # Sort and filter ranges
            min_start = wind_df['start'].min()
            max_end = wind_df['end'].max()

            extended_df = df[(df['start'] >= min_start) & (df['end'] <= max_end)].copy()
            start_idx = df[df['start'] == min_start].index[0]
            end_idx = df[df['end'] == max_end].index[0]

            start_extension = df.iloc[max(0, start_idx - n_wind):start_idx]
            end_extension = df.iloc[end_idx + 1:end_idx + 1 + n_wind]

            extended_wind_df = pd.concat([start_extension, extended_df, end_extension]).drop_duplicates().reset_index(drop=True)

            output_file_2 = os.path.join(tmp_folder, f'{pheno}_{ancestry}_{chr}_2.tsv')
            print(f"Writing header to {output_file_2}")
            awk_header_command = f"awk 'NR == 1 || NR == 2 {{print}}' {msp_file} > {output_file_2}"
            subprocess.run(awk_header_command, shell=True, check=True)

            for _, row in extended_wind_df.iterrows():
                awk_command = f"awk '$1 == \"{row['chr']}\" && $2 == {row['start']} && $3 == {row['end']} {{print}}' {msp_file} >> {output_file_2}"
                subprocess.run(awk_command, shell=True, check=True)

            
            laiobj = su.MSPReader(output_file_2).read()
            df_counts = pd.DataFrame(columns=['#IID', 'AFR', 'EAS', 'EUR'])
            # Initialize df_counts if empty
            if df_counts.empty:
                df_counts['#IID'] = laiobj.samples
                df_counts.iloc[:, 1:] = 0  # Initialize counts to 0

            # Update the DataFrame with counts
            n_wind, n_haplotypes = laiobj.lai.shape
            n_samp = len(laiobj.samples)
            print(laiobj.ancestry_map)

            print('start counting')

            # Vettorializzazione per evitare il ciclo annidato
            hap_to_sample = np.arange(n_haplotypes) // 2  # Map haplotype index to sample index
            ancestry_codes = laiobj.lai  # Entire LAI matrix


            ancestry_to_col = {v: i for i, v in enumerate(df_counts.columns[1:])}  # Map ancestry labels to columns

            temp_counts = np.zeros((n_samp, len(df_counts.columns) - 1), dtype=int)  # Temporary counts array

                # Vettorializzare il calcolo dei counts
            for wind in range(n_wind):
                ancestry_labels = [laiobj.ancestry_map.get(str(code), None) for code in ancestry_codes[wind, :]]
                ancestry_indices = [ancestry_to_col.get(label, None) for label in ancestry_labels]

                valid_indices = [(hap_idx, col) for hap_idx, col in enumerate(ancestry_indices) if col is not None]
                for hap_idx, col in valid_indices:
                    temp_counts[hap_to_sample[hap_idx], col] += 1

            # Aggiungiamo i counts nel df
            for i, sample in enumerate(laiobj.samples):
                df_counts.loc[i, df_counts.columns[1:]] = temp_counts[i]
            
            covar_df = pd.DataFrame(columns=['IID', f'{ancestry}', f'{ancestry}_percent'])

            # Calcolo dei total counts (uguale per tutti i campioni)
            total_counts = 2 * n_wind

            # Copia degli IID dal vecchio covar
            covar_df['IID'] = laiobj.samples

            # Estrazione dei counts per l'ancestry corrente
            ancestry_counts = df_counts[ancestry]

            # Calcolo della percentuale di ancestry
            covar_df[f'{ancestry}_percent'] = ancestry_counts / total_counts * 100

            # Ottieni il massimo dei counts per ogni sample (per tutti gli ancestries)
            max_counts_per_sample = df_counts.iloc[:, 1:].max(axis=1)

            # Aggiungi solo le due colonne necessarie e rimuovi la colonna binary
            covar_df[f'{ancestry}'] = (ancestry_counts == max_counts_per_sample).astype(int)

            # subset old covar
            old_covar_df['IID'] = old_covar_df['IID'].astype(str)
            covar_df['IID'] = covar_df['IID'].astype(str)

            old_covar_df = old_covar_df[old_covar_df['IID'].isin(covar_df['IID'])]

            # Unisci il vecchio covar_df con il nuovo (aggiungi solo le nuove colonne)
            covar_df = pd.merge(old_covar_df, covar_df[['IID', f'{ancestry}', f'{ancestry}_percent']], on='IID')

            covar_df.to_csv(os.path.join(output_folder, f'{ancestry}_{pheno}_covar_chr{chr}.tsv'), sep='\t', index=False)