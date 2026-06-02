library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

ANNOT_FILE  <- '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_analysis/credible_set_annotations.tsv'
OUTPUT_DIR  <- '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/gpn_msa_results/GO'
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ── Load genes ────────────────────────────────────────────────────────────────
annot <- read.table(ANNOT_FILE, sep='\t', header=TRUE, stringsAsFactors=FALSE)
genes <- unique(annot$gene_symbol)
genes <- genes[genes != '' & !is.na(genes)]
cat(sprintf('Gene symbols: %d\n', length(genes)))

# Convert symbols to Entrez IDs
entrez <- bitr(genes, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
cat(sprintf('Mapped to Entrez: %d/%d\n', nrow(entrez), length(genes)))

# ── GO enrichment (BP, MF, CC) ────────────────────────────────────────────────
ontologies <- c('BP', 'MF', 'CC')

for (ont in ontologies) {
    cat(sprintf('\nRunning GO %s...\n', ont))
    res <- enrichGO(
        gene          = entrez$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        pAdjustMethod = 'BH',
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.20,
        readable      = TRUE
    )

    if (is.null(res) || nrow(res) == 0) {
        cat(sprintf('  No significant terms for %s\n', ont))
        next
    }

    cat(sprintf('  Significant terms: %d\n', nrow(res)))

    # Save TSV
    write.table(as.data.frame(res),
                file=file.path(OUTPUT_DIR, sprintf('GO_%s_results.tsv', ont)),
                sep='\t', row.names=FALSE, quote=FALSE)

    # Dotplot (top 20 terms)
    p <- dotplot(res, showCategory=20) +
        ggtitle(sprintf('GO %s — SuSiE credible set genes', ont)) +
        theme(plot.title = element_text(size=12),
              axis.text.y = element_text(size=9))

    ggsave(file.path(OUTPUT_DIR, sprintf('GO_%s_dotplot.png', ont)),
           p, width=10, height=8, dpi=150)
    cat(sprintf('  Saved → GO_%s_dotplot.png\n', ont))
}

cat('\nDone.\n')
