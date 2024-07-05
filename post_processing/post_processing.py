from post_processing_functions import *
import sys

if __name__ == '__main__':
    # Check if the dataset argument is provided
    if len(sys.argv) != 2:
        print("Usage: python post_processing.py <dataset>")
        sys.exit(1)

    # The dataset variable is taken from the command line argument
    dataset = sys.argv[1]

    ancestry_list = ["AFR", "AHG", "EAS", "EUR", "NAT", "SAS", "WAS"]
    #ancestry_list = ["AFR"]

    # Use the dataset variable to construct file paths
    phe_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/{dataset}'
    general_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/{dataset}/output_ancestry_#/*/output.*.glm.logistic.hybrid'
    output_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/{dataset}/P_info_#.tsv'
    window_input_file = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/{dataset}/ancestry_AFR.vcf'
    wind_output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/{dataset}'
    plot_output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/plots/{dataset}'
    general_output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/{dataset}'
    FUMA_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/{dataset}'
    
    # extraction of the windows positions
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