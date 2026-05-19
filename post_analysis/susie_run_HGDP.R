library(susieR)

input_folder  <- '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_inputs_HGDP'
output_folder <- '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_results_HGDP'
plots_folder  <- '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_plots_HGDP'

dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)
dir.create(plots_folder,  recursive=TRUE, showWarnings=FALSE)

hits <- sort(list.dirs(input_folder, full.names=FALSE, recursive=FALSE))

all_results <- list()

for (hit in hits) {
    cat(sprintf("\n%s\n", hit))

    z_file  <- file.path(input_folder, hit, paste0(hit, '_zscores.tsv'))
    ld_file <- file.path(input_folder, hit, paste0(hit, '_ld.gz'))

    if (!file.exists(z_file) || !file.exists(ld_file)) {
        cat("  Missing files, skipping\n")
        next
    }

    zscores <- read.table(z_file, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    z <- zscores$Z
    n <- as.integer(median(zscores$N))
    p <- length(z)

    cat(sprintf("  SNPs: %d  N: %d\n", p, n))

    cat("  Loading LD matrix...\n")
    R <- as.matrix(read.table(gzfile(ld_file), header=FALSE))
    R <- (R + t(R)) / 2

    if (nrow(R) != p) {
        cat(sprintf("  Dimension mismatch: LD %dx%d vs z %d, skipping\n", nrow(R), ncol(R), p))
        next
    }

    if (n < p) {
        cat(sprintf("  WARNING: N=%d < p=%d — LD rank-deficient, results may be unreliable\n", n, p))
    }

    cat("  Running SuSiE...\n")
    result <- tryCatch(
        susie_rss(z=z, R=R, n=n, L=10, verbose=FALSE,
                  estimate_residual_variance=FALSE),
        error = function(e) {
            cat(sprintf("  SuSiE error: %s\n", e$message))
            NULL
        }
    )

    if (is.null(result)) next

    pip_df <- data.frame(
        ID    = zscores$ID,
        CHROM = zscores$CHROM,
        POS   = zscores$POS,
        Z     = zscores$Z,
        PIP   = result$pip,
        stringsAsFactors = FALSE
    )
    pip_df <- pip_df[order(-pip_df$PIP), ]
    pip_nonaffx <- pip_df[!startsWith(pip_df$ID, 'Affx'), ]

    cs_list <- result$sets$cs
    cs_df   <- data.frame()
    n_cs    <- 0

    if (!is.null(cs_list) && length(cs_list) > 0) {
        n_cs <- length(cs_list)
        for (i in seq_along(cs_list)) {
            idx <- cs_list[[i]]
            cs_df <- rbind(cs_df, data.frame(
                CS    = i,
                ID    = zscores$ID[idx],
                CHROM = zscores$CHROM[idx],
                POS   = zscores$POS[idx],
                Z     = zscores$Z[idx],
                PIP   = result$pip[idx],
                stringsAsFactors = FALSE
            ))
        }
        top_row <- if (nrow(pip_nonaffx) > 0) pip_nonaffx[1, ] else pip_df[1, ]
        cat(sprintf("  Credible sets: %d\n", n_cs))
        cat(sprintf("  Top PIP: %.4f  (%s  chr%s:%s)\n",
                    top_row$PIP, top_row$ID, top_row$CHROM, top_row$POS))
    } else {
        cat("  No credible sets found\n")
    }

    hit_out <- file.path(output_folder, hit)
    dir.create(hit_out, recursive=TRUE, showWarnings=FALSE)

    write.table(pip_df,
                file=file.path(hit_out, paste0(hit, '_pip.tsv')),
                sep='\t', row.names=FALSE, quote=FALSE)

    if (nrow(cs_df) > 0) {
        write.table(cs_df,
                    file=file.path(hit_out, paste0(hit, '_credible_sets.tsv')),
                    sep='\t', row.names=FALSE, quote=FALSE)
    }

    top_row <- if (nrow(pip_nonaffx) > 0) pip_nonaffx[1, ] else pip_df[1, ]
    all_results[[hit]] <- data.frame(
        hit     = hit,
        n_snps  = p,
        n       = n,
        n_cs    = n_cs,
        top_snp = if (n_cs > 0) top_row$ID  else NA,
        top_pip = if (n_cs > 0) top_row$PIP else NA,
        stringsAsFactors = FALSE
    )

    png(file.path(plots_folder, paste0(hit, '_pip.png')), width=1400, height=500)
    susie_plot(result, y='PIP',
               main=sprintf('%s  |  CS: %d', hit, n_cs))
    dev.off()

    png(file.path(plots_folder, paste0(hit, '_z.png')), width=1400, height=500)
    susie_plot(result, y='z',
               main=sprintf('%s  |  z-scores', hit))
    dev.off()

    cat(sprintf("  Saved → %s/\n", hit_out))
}

cat(sprintf("\n%s\n", strrep("=", 60)))
if (length(all_results) > 0) {
    summary_df   <- do.call(rbind, all_results)
    summary_file <- file.path(output_folder, 'susie_summary_HGDP.tsv')
    write.table(summary_df, file=summary_file, sep='\t', row.names=FALSE, quote=FALSE)
    cat(sprintf("Summary saved → %s\n", summary_file))
    print(summary_df)
}
cat("Done.\n")
