#!/usr/bin/env -S Rscript --vanilla

# anndata_concat.R
# This script concatenates multiple AnnData files into a single h5ad object.
#
# Usage:
#   Rscript anndata_concat.R -i <input> -o <output_file>

suppressPackageStartupMessages(library(optparse))

# Define command-line options for the script
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "Input h5ad file",
              metavar = "character"),
  make_option(c("-o", "--output_file"),
              type = "character",
              default = NULL,
              help = "Output coloc guide table",
              metavar = "character")
)
# Parse the command-line options
opt_parser <- OptionParser(usage = "Usage: %prog -i <input> -o <output_file>",
                           option_list = option_list)
opt <- parse_args(opt_parser)

# Read list of .h5ad files
sce <- zellkonverter::readH5AD(opt$input,reader="R")
coloc_input <- flanders::anndata2coloc_input(sce)
data.table::fwrite(coloc_input, file = opt$output_file)