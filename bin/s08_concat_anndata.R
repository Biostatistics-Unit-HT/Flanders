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
              help = "File containin a list of .h5ad files",
              metavar = "character"),
  make_option(c("-o", "--output_file"),
              type = "character",
              default = NULL,
              help = "Name of concatenated .h5ad output file",
              metavar = "character")
)
# Parse the command-line options
opt_parser <- OptionParser(usage = "Usage: %prog -i <input> -o <output_file>",
                           option_list = option_list)
opt <- parse_args(opt_parser)

# Read list of .h5ad files
input_files <- readLines(opt$input)
list_of_ads <- lapply(input_files,anndata::read_h5ad)

# Concatenate the AnnData objects
merged_ad <- anndata::concat(
	list_of_ads,
	join="outer",
	merge="first"
)
merged_ad <- flanders::fix_ad_var(merged_ad)

# Write the merged AnnData object to a .h5ad file
anndata::write_h5ad(merged_ad, filename = opt$output_file)