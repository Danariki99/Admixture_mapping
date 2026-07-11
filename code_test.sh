
#!/bin/bash
# =============================================================================
# END-TO-END TEST PIPELINE  —  NEW pipeline
# =============================================================================
# Walks the whole current pipeline on a small test dataset, coherent with the
# production pipeline we now use:
#
#   LAI (GNomix, unchanged)  ->  admixture mapping (NEW covariates + NEW
#   significance criteria)  ->  CONDITIONAL fine mapping (SNP ADD + LAI in the
#   same model)  ->  Borzoi SED on ALL window SNPs  ->  example best-SED plot.
#
# NOTES
#  * LAI is still GNomix (the RFMix swap is planned separately). All GNomix usage
#    is left exactly as it is and can be used as-is.
#  * Heavy / cluster / GPU steps are shown but COMMENTED (like the GNomix step):
#    they submit SLURM jobs or need the separate Borzoi environment. Uncomment /
#    adapt them for a real run.
#  * Steps whose *_test version does not exist yet reference the PRODUCTION
#    script (hardcoded paths). Either add a *_test variant or point its paths at
#    the test dataset — left for us to finish where needed.
# =============================================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <data_folder> <result_folder>"
    exit 1
fi

data_folder="$1"
result_folder="$2"

vcf_folder="${result_folder}/vcf_folder"
msp_folder="${result_folder}/msp_folder"

# =============================================================================
# 1) LOCAL ANCESTRY INFERENCE  (GNomix — unchanged)
# =============================================================================
#python LAI/chrom_division.py --results "$result_folder" --data "$data_folder"
echo "Using GnoMix software for LAI"
#python LAI/gnomix_training_test.py --vcf_folder "$vcf_folder" --data_folder "$data_folder" --result_folder "$result_folder"

./LAI/files_moving.sh "$result_folder"

# =============================================================================
# 2) ADMIXTURE MAPPING
# =============================================================================
python pre_processing/pre_processing_test.py "$msp_folder" "$result_folder"

# local-ancestry logistic association, per ancestry / phenotype
./association_execution/analysis_execution_test.sh "$result_folder" "$data_folder"

# post-processing with the NEW significance criteria:
#   BY (fdr_by) correction, >=5 contiguous windows, <=1 expected false positive
python post_processing/post_processing_test.py "$result_folder" "$data_folder"

python post_analysis/ancestry_counts_test.py "$msp_folder" "$result_folder"

# =============================================================================
# 3) COVARIATES  (NEW: global-ancestry PROPORTIONS + PCs)
# =============================================================================
python post_analysis/keep_files_test.py "$result_folder"

python post_analysis/keep_files_post_processing_test.py "$result_folder" "$data_folder"

python post_analysis/LD_pruning_test.py "$result_folder" "$data_folder"

python post_analysis/PCA_generation_test.py "$result_folder" "$data_folder"

# NEW covariates: global-ancestry proportions (not the old hard-call covariates)
# + PCs. If the test covar files are not regenerated here, keep the existing file
# and edit it (the production version is post_analysis/new_covariable_creation.py).
python post_analysis/new_covariable_creation_test.py "$result_folder" "$data_folder"

python post_analysis/PCA_covar_files_creation_test.py "$result_folder"

# =============================================================================
# 4) FINE MAPPING  (NEW: CONDITIONAL analysis — SNP ADD + LAI covariate)
# =============================================================================
python post_vcf/keep_creation_chrom_test.py "$result_folder" "$data_folder"

python post_vcf/snps_creation_more_wind_test.py "$result_folder" "$data_folder"

# per-window covariate files that INCLUDE the window's LAI (local ancestry) as a
# covariate, so the GLM tests the SNP (ADD) and the ancestry (LAI) jointly.
# Covar columns: IID age sex BMI AFR AHG EAS EUR NAT OCE SAS WAS LAI
python post_analysis/window_covar_creation_test.py "$result_folder" "$data_folder"

# conditional fine mapping: plink2 --glm firth-fallback with --covar-col-nums
# 2-7,9-13 (EUR = col 8 dropped as reference, LAI = col 13) -> ADD + LAI test rows
./post_analysis/fine_mapping_new_test.sh "$result_folder" "$data_folder"

# fine-mapping post-processing: BY correction on ADD and LAI p-values;
# candidates = ADD significant AND LAI NOT significant AND OR direction concordant
# with the admixture-mapping OR; then extract the causal candidate per hit.
python post_analysis/fine_mapping_post_processing_test.py "$result_folder" "$data_folder"
python post_analysis/fine_mapping_causal_candidate_extraction_test.py "$result_folder"

# =============================================================================
# 5) BORZOI  (SED on ALL window SNPs)  +  example best-SED plot
# =============================================================================
# build the hg38 VCF of EVERY SNP tested in the fine-mapping windows (lifts hg19->hg38)
python post_analysis/windows_snps_vcf_creation_test.py "$result_folder"

# run Borzoi SED (gene-level) across the 4 model folds. Needs the SEPARATE borzoi
# environment (TensorFlow) and, in practice, a GPU. See singularity.def — kept
# commented so a plain (CPU / no-Borzoi) test run does not fail here.
#python post_analysis/borzoi_sed_windows_test.py "$result_folder"

# average the 4 folds -> per-(SNP, gene) SED tables
#python post_analysis/borzoi_sed_post_processing_test.py "$result_folder"

# EXAMPLE plot of the best SED results (top |SED| SNP-gene pairs) — illustrative
# only, not the full production figure set.
python post_analysis/borzoi_sed_example_plot.py "$result_folder"
