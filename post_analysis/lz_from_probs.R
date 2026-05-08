
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(locuszoomr)
  library(htmlwidgets)
  library(plotly)
})

opt_list <- list(
  make_option("--metal", type="character", help="TSV with CHR,BP,SNP,P"),
  make_option("--refsnp", type="character", help="Lead SNP (rsID or CHR:POS)"),
  make_option("--build", type="character", default="GRCh37"),
  make_option("--pop", type="character", default="EUR",
              help="1000G population for LD: EUR/AFR/EAS/SAS/AMR/ALL"),
  make_option("--out_html", type="character", help="Output HTML file"),
  make_option("--flank", type="integer", default=400000, help="+/- bp window"),
  make_option("--ldlink_token", type="character", default=Sys.getenv("LDLINK_TOKEN"))
)
opt <- parse_args(OptionParser(option_list=opt_list))
if (is.null(opt$metal) || is.null(opt$refsnp) || is.null(opt$out_html)) {
  stop("Missing required args: --metal, --refsnp, --out_html")
}
if (opt$ldlink_token == "") {
  stop("LDLINK_TOKEN not set. export LDLINK_TOKEN=... before running.")
}

# Read association table (must contain CHR,BP,SNP,P)
assoc <- fread(opt$metal, sep="\t", header=TRUE, showProgress=FALSE)
setnames(assoc, old=c("CHR","BP","SNP","P"),
                new=c("CHR","BP","SNP","P"), skip_absent=TRUE)
assoc[, CHR := as.character(CHR)]
assoc[, BP  := as.integer(BP)]

# Build locus object centered on refsnp
L <- locus(assoc = assoc, build = opt$build,
           refsnp = opt$refsnp, flank = opt$flank)

# Fetch LD from LDlink (1000G) via API token
# e.g. pop: EUR/AFR/EAS/SAS/AMR/ALL
Sys.setenv(LDLINK_TOKEN = opt$ldlink_token)
L <- link_LD(L, population = opt$pop)  # adds r2 for coloring

# Interactive plot with genes + recombination
p <- locus_plotly(L, show_genes = TRUE, show_recomb = TRUE)

# Write HTML (self-contained)
dir.create(dirname(opt$out_html), recursive = TRUE, showWarnings = FALSE)
saveWidget(as_widget(p), opt$out_html, selfcontained = TRUE)
cat("Wrote", opt$out_html, "\n")
