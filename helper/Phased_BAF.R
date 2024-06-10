#!/usr/bin/env Rscript

# Load necessary libraries
libraries <- c("optparse", "dplyr", "VariantAnnotation", "GenomicRanges", "data.table", "pbmcapply")
for (lib in libraries) {
  suppressMessages({
    if (!require(lib, character.only = TRUE)) {
      stop(paste("The package", lib, "is not installed."))
    }
  })
}


source("mod/mod_allele_imbalance.R")

# Define command line options
option_list <- list(
  make_option(c("-I", "--input"), type = "character", default = NULL, help = "Path to the joint genotyped VCF file", metavar = "character"),
  make_option(c("-O", "--output"), type = "character", default = NULL, help = "Path to the output file, if NULL stdout", metavar = "character"),
  make_option(c("-C", "--chr"), type = "character", default = NULL, help = "chr region for analysis, e.g. chr1, all, auto", metavar = "character"),
  make_option(c("-R", "--ref"), type = "character", default = "hg38", help = "Reference genome, e.g. hg38, hg19, [default -R hg38]", metavar = "character"),
  make_option(c("-N", "--cores"), type = "character", default = "1", help = "Number of cores for parallel analysis, [default -N 1]", metavar = "character")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if no options provided and print help
if (is.null(opt$input) && is.null(opt$output) && is.null(opt$chr) && opt$ref == "hg38" && opt$cores == "1") {
  print_help(opt_parser)
  quit(status = 0)  # Exit cleanly if no options are provided
}

# Validate essential inputs
if (is.null(opt$input) || is.null(opt$chr) ) {
  print_help(opt_parser)
  stop("Both --input, --output and --chr must be provided")
}

# Ensure input file exists
if (!file.exists(opt$input)) {
  stop("Input file does not exist.")
}

# Determine the output file path
output_path <- ifelse(grepl("^/", opt$output),
                      paste0(opt$output, ".phased_snp_summary.tsv"),  # Absolute path
                      paste0(getwd(), "/", opt$output, ".phased_snp_summary.tsv"))  # Relative path

# Handle special chromosome specifications
if (opt$chr %in% c("all", "auto")) {
  chr_map <- list(all = paste0("chr", c(1:22, "X")), auto = paste0("chr", 1:22))
  opt$chr <- chr_map[[opt$chr]]
} else {
  opt$chr <- opt$chr
}

# Load reference genome information
ref_file_path <- sprintf("data/%s.info.txt", opt$ref)
if (!file.exists(ref_file_path)) {
  stop("Reference file does not exist.")
}
ref_info <- fread(ref_file_path)


tryCatch({
  # Parallel processing of SNPs
  result <- pbmclapply(opt$chr, function(chr) {
    loc.end <- ref_info %>% filter(chrom == chr | chrom == paste0("chr", chr)) %>% pull(seqlengths)
    range.gr <- GRanges(chr, ranges = IRanges(start = 1, end = loc.end))
    snp.gr <- ReadGVCF(opt$input, opt$ref, range.gr, target_spacing = 0)
    snp.df <- as.data.frame(snp.gr)
    snp.df %>% group_by(B_InhFrom) %>%
      summarise(count = n(), .groups = 'drop') %>%
      mutate(total = sum(count), freq = round(count / total,digits = 3), chrom = chr)
  }, mc.cores = as.numeric(opt$cores))
  
  snp_stats <- rbindlist(result)
  rm(result); gc()
  # Conditional output based on whether -O is provided
  if (is.null(opt$output)) {
    # No output file specified, write to stdout
    write.table(snp_stats, file="", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  } else {
    # Output file specified, write to file
    output_path <- ifelse(grepl("^/", opt$output),
                          opt$output,  # Absolute path provided
                          paste0(getwd(), "/", opt$output))  # Relative path provided
    write.table(snp_stats, file=output_path, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  }
}, error = function(e) {
  cat("An error occurred: ", e$message, "\n")
})