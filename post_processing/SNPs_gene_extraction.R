#!/usr/bin/env Rscript

# deps
require(biomaRt)
require(data.table)
require(optparse)

# Define command line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the file with SNP positions", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

# Parse command line arguments
parser = OptionParser(option_list=option_list)
args = parse_args(parser)

# Ensure the output directory ends with a slash
output_dir <- normalizePath(args$output, mustWork = FALSE, winslash = "/")
if (!grepl("/$", output_dir)) {
  output_dir <- paste0(output_dir, "/")
}

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
}

# Use the input file path directly
wind.file <- args$input

# Ensure the file exists
if (!file.exists(wind.file)) {
  stop("Input file does not exist: ", wind.file)
}

# Read important positions
start <- fread(wind.file)$POS
end <- fread(wind.file)$end_POS
chromosome <- fread(wind.file)$`#CHROM`

# Connect to dataset
mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "https://feb2014.archive.ensembl.org")

# Empty data frame to store ensembl outputs
snp.tab <- data.frame()

for (i in 1:length(start)) {
  print(paste("Retrieving SNPs for position ", start[i], " to ", end[i], " on chromosome ", chromosome[i], " (", i, "/", length(start), ")", sep = ""))

  # Create coords vector as required for getBM
  query <- paste(chromosome[i], start[i], end[i], sep = ":")

  # set maximum and current number of attempts
  max_attempts <- 5
  current_attempt <- 1

  # try X times per SNP set
  while (current_attempt <= max_attempts) {
    tryCatch({
      print(paste("Attempt ", current_attempt, " of gene ", i, "...", sep = ""))

      # Get table from ensembl
      sub.snp.tab <- getBM(
        attributes = c("chr_name", "chrom_start", "refsnp_id", "allele"),
        filters = c("chromosomal_region"), values = query, mart = mart
      )

      # If the command succeeds, break out of the loop
      break
    }, error = function(err) {
      # Print the error message
      cat(paste("Attempt", current_attempt, "failed with error:", conditionMessage(err), "\n"))

      # Increment the attempt counter
      current_attempt <- current_attempt + 1

      # If it's the last attempt, stop trying
      if (current_attempt > max_attempts) {
        stop("Max attempts reached. Exiting.")
      }
    })
  }

  # skip if number of SNPs is 0
  if (nrow(sub.snp.tab) == 0) next

  # Append to general snp tab
  snp.tab <- rbind(snp.tab, sub.snp.tab)
}

# Store table
output_file <- file.path(output_dir, "significant_SNPs.txt")
write.table(snp.tab, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE)

# Print the output file path
cat(output_file)