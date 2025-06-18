from post_processing_functions_test import *
import os
import sys

if __name__ == '__main__':
    result_folder = sys.argv[1]
    data_folder = sys.argv[2]
    # Lista ancestry usate nei file
    ancestry_list = ["AFR", "EAS", "EUR"]

    # Nuovi percorsi coerenti con il setup "data/" e "results/"
    phe_folder = os.path.join(data_folder, 'phe_files')
    window_input_file = os.path.join(result_folder, 'vcf_files', 'ancestry_AFR.vcf')
    wind_output_folder = os.path.join(result_folder, 'post_processing_files')
    plot_output_folder = os.path.join(result_folder, 'plots')
    general_output_folder = os.path.join(result_folder, 'post_processing_files')
    FUMA_folder = os.path.join(result_folder, 'FUMA')    

    os.makedirs(wind_output_folder, exist_ok=True)
    os.makedirs(plot_output_folder, exist_ok=True)
    os.makedirs(general_output_folder, exist_ok=True)
    os.makedirs(FUMA_folder, exist_ok=True)

    # Estrazione delle posizioni delle finestre
    print("Extracting windows positions")
    window_pos_file = positions_extraction(window_input_file, wind_output_folder)

    for ancestry in ancestry_list:
        print(f"\nProcessing ancestry: {ancestry}")

        # Percorso generico degli output PLINK per questa ancestry
        general_file = os.path.join(result_folder, f'output_ancestry_{ancestry}/*/output.*.glm.logistic.hybrid')
        output_file = os.path.join(result_folder, f'post_processing_files/P_info_{ancestry}.tsv')

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
