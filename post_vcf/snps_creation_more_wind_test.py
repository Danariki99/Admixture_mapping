import os
import sys
import subprocess
import shutil
import pandas as pd

if __name__ == '__main__':
    # Percorsi coerenti con la pipeline test
    wind_folder = './results/FUMA'
    vcf_folder = './data/vcf_files'
    output_folder = './results/snps_files'
    tmp_folder = './results/tmp'

    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(tmp_folder, exist_ok=True)

    list_of_files = [f for f in os.listdir(wind_folder) if f.endswith(".txt")]

    for wind_filename in list_of_files:
        ancestry = wind_filename.split('_')[0]
        input_file = os.path.join(vcf_folder, f"chr1_{ancestry}.vcf")  # usa chr1 come riferimento per posizione
        wind_file = os.path.join(wind_folder, wind_filename)
        wind_df_1 = pd.read_csv(wind_file, sep='\t')

        chroms = wind_df_1['chr'].unique()

        for chr in chroms:
            output_file_1 = os.path.join(tmp_folder, 'chrom_pos_tmp.tsv')

            # estrai CHROM e POS dalla VCF ancestry-specifico
            awk_command = f"awk '!/^##/ && $1 == \"{chr}\" {{print $1, $2}}' {input_file} > {output_file_1}"
            subprocess.run(awk_command, shell=True, check=True)

            df = pd.read_csv(output_file_1, sep=' ', header=None, names=['chr', 'pos'])
            df[['start', 'end']] = df['pos'].astype(str).str.split('_', expand=True)
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            df = df.drop(columns=['pos'])

            wind_df = wind_df_1.loc[wind_df_1['chr'] == chr]

            n_wind = 3  # estensione prima e dopo
            df = df.sort_values(['start']).reset_index(drop=True)
            wind_df = wind_df.sort_values(['start']).reset_index(drop=True)

            min_start = wind_df['start'].min()
            max_end = wind_df['end'].max()

            extended_df = df[(df['start'] >= min_start) & (df['end'] <= max_end)].copy()
            start_idx = df[df['start'] == min_start].index[0]
            end_idx = df[df['end'] == max_end].index[0]

            start_extension = df.iloc[max(0, start_idx - n_wind):start_idx]
            end_extension = df.iloc[end_idx + 1:end_idx + 1 + n_wind]

            extended_wind_df = pd.concat([start_extension, extended_df, end_extension]).drop_duplicates().reset_index(drop=True)

            output_file = os.path.join(output_folder, wind_filename.replace('wind', f'snps_chr{chr}'))
            tmp_file_vcf = os.path.join(tmp_folder, wind_filename.replace('wind.txt', 'tmp.vcf'))
            tmp_file_snps = os.path.join(tmp_folder, wind_filename.replace('wind', 'tmp_snps'))

            list_of_snps = []

            for _, row in extended_wind_df.iterrows():
                chrom = row['chr']
                start = row['start']
                end = row['end']

                plink_command = [
                    "/private/home/rsmerigl/plink2",
                    "--vcf", input_file,
                    "--chr", str(chrom),
                    "--from-bp", str(start),
                    "--to-bp", str(end),
                    "--recode", "vcf",
                    "--out", tmp_file_vcf.replace('.vcf', '')
                ]
                subprocess.run(plink_command, check=True)

                awk_command = f"awk '!/^##/ {{print $3}}' {tmp_file_vcf} > {tmp_file_snps}"
                subprocess.run(awk_command, shell=True, check=True)

                df_snps = pd.read_csv(tmp_file_snps, header=None, names=['ID'])
                list_of_snps.extend(df_snps['ID'].dropna().tolist())

                # reset temporanei
                shutil.rmtree(tmp_folder)
                os.makedirs(tmp_folder)

            last_df = pd.DataFrame(list(set(list_of_snps)), columns=['#ID'])
            last_df.to_csv(output_file, index=False)
            print(f"✅ File {output_file} (chr {chrom}) creato")
