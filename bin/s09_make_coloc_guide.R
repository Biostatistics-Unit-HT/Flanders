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
  make_option(c("-x", "--exclude"),
              type = "character",
              default = NULL,
              help = "Tab separated file with header study_id, phenotype_id that report traits to exclude",
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

# Read h5ad file and make guide table
message("Reading h5ad file")
sce <- zellkonverter::readH5AD(opt$input,reader="R")

message("Making coloc guide table")
coloc_input <- flanders::anndata2coloc_input(sce)
message(nrow(coloc_input), " coloc tests generated")

# Read the exclusion table using read_delim and convert to data.table
exclude_dt <- readr::read_delim(opt$exclude, delim = "\t", col_names = TRUE)
exclude_dt <- as.data.table(exclude_dt)
exclude_dt <- exclude_dt[study_id != "EMPTY"]

# If there are studies to exclude, remove them
if (nrow(exclude_dt) > 0) {
  message("Found ", nrow(exclude_dt), " studies to exclude")
  exclude_dt[, unique_id := paste(study_id, phenotype_id, sep = "_")]
  coloc_input[, t1_uniqueid := paste(t1_study_id, t1_phenotype_id, sep = "_")]
  coloc_input[, t2_uniqueid := paste(t2_study_id, t2_phenotype_id, sep = "_")]
  coloc_input[, t1_excluded := t1_uniqueid %in% exclude_dt$unique_id]
  coloc_input[, t2_excluded := t2_uniqueid %in% exclude_dt$unique_id]
  coloc_input <- coloc_input[!(coloc_input$t1_excluded & coloc_input$t2_excluded)]

  message(nrow(coloc_input), " coloc tests remaining after exclusion")
}

data.table::fwrite(coloc_input, file = opt$output_file)
message("Coloc guide table saved to ", opt$output_file)
