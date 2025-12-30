import argparse
from post_processing_functions import *

def parse_args():
    parser = argparse.ArgumentParser(description="Run admixture-mapping post processing steps.")
    parser.add_argument("dataset", help="Dataset identifier (e.g. ukbb, all_of_us)")
    parser.add_argument(
        "--apply-fb-qc",
        action="store_true",
        help="Drop windows whose RFMix FB max-confidence is below the threshold."
    )
    parser.add_argument(
        "--fb-template",
        default=None,
        help="Optional template for FB files. Use {dataset}, {ancestry} and {chrom} placeholders."
    )
    parser.add_argument(
        "--fb-root",
        default=None,
        help="Fallback directory where FB files are stored when --fb-template is not provided."
    )
    parser.add_argument(
        "--fb-threshold",
        type=float,
        default=0.9,
        help="Minimum max-confidence required to keep a window (default: 0.9)."
    )
    parser.add_argument(
        "--fb-chunksize",
        type=int,
        default=200000,
        help="Number of rows read per chunk when scanning FB files (default: 200000)."
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    dataset = args.dataset

    ancestry_list = ["AFR", "AHG", "EAS", "EUR", "NAT", "SAS", "WAS"]

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
    for ancestry in ancestry_list:
    # result analysis and manhattan plot creation
        print("Analyzing results")
        significant_position_file = result_analysis(
            [ancestry],
            phe_folder,
            general_file,
            window_pos_file,
            output_file,
            plot_output_folder,
            general_output_folder,
            dataset_name=dataset,
            apply_fb_filter=args.apply_fb_qc,
            fb_template=args.fb_template,
            fb_root=args.fb_root,
            fb_threshold=args.fb_threshold,
            fb_chunksize=args.fb_chunksize
        )
        if significant_position_file is not None:
            # extraction of the SNPs
            print("Extracting SNPs")
            total_SNPs_file = SNPs_extraction(significant_position_file, general_output_folder)

            # associate SNPs to windows
            print("Associating SNPs to windows")
            SNPs_with_P_info = associate_SNPs_to_windows(total_SNPs_file, significant_position_file, general_output_folder)

            # create FUMA ready files
            print("Creating FUMA ready files")
            output_folder_snps, output_folder_wind = FUMA_files_creation(SNPs_with_P_info, FUMA_folder)
            
