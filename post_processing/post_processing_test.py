from post_processing_functions import *
import os

if __name__ == '__main__':
    # Lista ancestry usate nei file
    ancestry_list = ["AFR", "AHG", "EAS", "EUR", "NAT", "SAS", "WAS"]

    # Nuovi percorsi coerenti con il setup "data/" e "results/"
    phe_folder = './data/phe_files'
    window_input_file = './data/vcf_files/ancestry_AFR.vcf'  # puoi cambiarlo con qualsiasi .vcf valido
    wind_output_folder = './results/post_processing_files'
    plot_output_folder = './results/plots'
    general_output_folder = './results/post_processing_files'
    FUMA_folder = './results/FUMA'

    # Estrazione delle posizioni delle finestre
    print("Extracting windows positions")
    window_pos_file = positions_extraction(window_input_file, wind_output_folder)

    for ancestry in ancestry_list:
        print(f"\nProcessing ancestry: {ancestry}")

        # Percorso generico degli output PLINK per questa ancestry
        general_file = f'./results/output_ancestry_{ancestry}/*/output.*.glm.logistic.hybrid'
        output_file = f'./results/post_processing_files/P_info_{ancestry}.tsv'

        # Analisi dei risultati e creazione dei plot
        print("Analyzing results")
        significant_position_file = result_analysis(
            [ancestry],
            phe_folder,
            general_file,
            window_pos_file,
            output_file,
            plot_output_folder,
            general_output_folder
        )

        if significant_position_file is not None:
            print("Extracting SNPs")
            total_SNPs_file = SNPs_extraction(significant_position_file, general_output_folder)

            print("Associating SNPs to windows")
            SNPs_with_P_info = associate_SNPs_to_windows(total_SNPs_file, significant_position_file, general_output_folder)

            print("Creating FUMA ready files")
            output_folder_snps, output_folder_wind = FUMA_files_creation(SNPs_with_P_info, FUMA_folder)
