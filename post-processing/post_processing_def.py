from post_processing_functions import *

if __name__ == '__main__':

    ancestry_list = ["AFR"]

    phe_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/ukbb'
    general_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb/output_ancestry_#/*/output.*.glm.logistic.hybrid'
    output_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/ukbb/P_info_#.tsv'
    window_input_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/ukbb/ancestry_AFR.vcf'
    wind_output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/ukbb'
    plot_output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/plots/ukbb'
    general_output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/ukbb'
    FUMA_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb'
    
    # extraction of the windows positiojns
    print("Extracting windows positions")
    window_pos_file = positions_extraction(window_input_file, wind_output_folder)

    # result analysis and manhattan plot creation
    print("Analyzing results")
    significant_position_file = result_analysis(ancestry_list, phe_folder, general_file, window_pos_file, output_file, plot_output_folder, general_output_folder)
    
    # extraction of the SNPs
    print("Extracting SNPs")
    total_SNPs_file = SNPs_extraction(significant_position_file, general_output_folder)

    # associate SNPs to windows
    print("Associating SNPs to windows")
    SNPs_with_P_info = associate_SNPs_to_windows(total_SNPs_file, significant_position_file, general_output_folder)

    # create FUMA ready files
    print("Creating FUMA ready files")
    output_folder_snps, output_folder_wind = FUMA_files_creation(SNPs_with_P_info, FUMA_folder)
