#!/usr/bin/env Rscript

options(timeout = 300)  # 5-minute HTTP timeout for Ensembl queries

# deps
require(biomaRt)
require(data.table)
require(optparse)

message(paste("R version:", R.version$major, R.version$minor))
message(paste("biomaRt version:", packageVersion("biomaRt")))

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

# Connect to Ensembl with retry across hosts
connect_mart <- function() {
  hosts <- c("https://grch37.ensembl.org", "https://useast.ensembl.org", "https://asia.ensembl.org")
  max_conn_attempts <- 5
  for (conn_attempt in 1:max_conn_attempts) {
    for (h in hosts) {
      m <- tryCatch(
        useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = h),
        error = function(e) { message(paste("Host", h, "failed:", conditionMessage(e))); NULL }
      )
      if (!is.null(m)) { message(paste("Connected to:", h)); return(m) }
    }
    if (conn_attempt < max_conn_attempts) {
      message(paste("All hosts failed — waiting 60s before retry", conn_attempt, "of", max_conn_attempts))
      Sys.sleep(60)
    }
  }
  stop("All Ensembl hosts failed after all retries.")
}

mart <- connect_mart()

# Empty data frame to store ensembl outputs
snp.tab <- data.frame()
skipped_windows <- c()

for (i in 1:length(start)) {
  message(paste0("Retrieving SNPs for position ", start[i], " to ", end[i],
                 " on chromosome ", chromosome[i], " (", i, "/", length(start), ")"))

  query <- paste(chromosome[i], start[i], end[i], sep = ":")

  max_attempts <- 5
  current_attempt <- 1
  skip_window <- FALSE
  sub.snp.tab <- NULL

  while (current_attempt <= max_attempts) {
    tryCatch({
      sub.snp.tab <- getBM(
        attributes = c("chr_name", "chrom_start", "refsnp_id", "allele"),
        filters = "chromosomal_region", values = query, mart = mart
      )
      break
    }, error = function(err) {
      message(paste0("  Attempt ", current_attempt, " failed: ", conditionMessage(err)))
      current_attempt <<- current_attempt + 1
      if (current_attempt > max_attempts) {
        message(paste0("  Window ", i, " (", query, ") failed after ", max_attempts, " attempts — skipping"))
        skip_window <<- TRUE
      } else {
        Sys.sleep(min(10 * current_attempt, 60))  # progressive backoff up to 60s
      }
    })
    if (skip_window) break
  }

  if (skip_window || is.null(sub.snp.tab) || nrow(sub.snp.tab) == 0) {
    if (skip_window) skipped_windows <- c(skipped_windows, query)
    next
  }

  snp.tab <- rbind(snp.tab, sub.snp.tab)
}

if (length(skipped_windows) > 0) {
  message(paste0("WARNING: ", length(skipped_windows), " windows skipped due to Ensembl errors:"))
  for (w in skipped_windows) message(paste0("  ", w))
}

# Store table
output_file <- file.path(output_dir, "significant_SNPs.txt")
write.table(snp.tab, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE)

# Print the output file path
cat(output_file)
