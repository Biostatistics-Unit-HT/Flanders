#!/usr/bin/env -S Rscript --vanilla

# anndata_concat.R
# This script concatenates multiple AnnData files into a single h5ad object.
#
# Usage:
#   Rscript anndata_concat.R -i <input> -o <output_file>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SingleCellExperiment))

# Define command-line options for the script
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "Input h5ad file",
              metavar = "character"),
  make_option(c("-b", "--exclude_bothsides"),
              type = "character",
              default = NULL,
              help = "Tab separated file with header study_id, phenotype_id that report traits to exclude when both t1 and t2 are in the list",
              metavar = "character"),
  make_option(c("-x", "--exclude_oneside"),
              type = "character",
              default = NULL,
              help = "Tab separated file with header study_id, phenotype_id that report traits to exclude when either t1 or t2 are in the list",
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
coloc_input <- as.data.table(flanders::anndata2coloc_input(sce))
message(nrow(coloc_input), " coloc tests generated")

# Make a unique id useful for filtering
coloc_input[, t1_uniqueid := paste(t1_study_id, t1_phenotype_id, sep = "_")]
coloc_input[, t2_uniqueid := paste(t2_study_id, t2_phenotype_id, sep = "_")]

if (!is.null(opt$exclude_bothsides)) {
  exclude_2sides <- readr::read_delim(opt$exclude_bothsides, delim = "\t", col_names = TRUE)
  exclude_2sides <- as.data.table(exclude_2sides)

  message("Found ", nrow(exclude_2sides), " entries in both side exclusion")
  exclude_2sides[, unique_id := paste(study_id, phenotype_id, sep = "_")]

  # Exclude studies where both t1 and t2 are in the exclusion list
  coloc_input[, t1_excluded := t1_uniqueid %in% exclude_dt$unique_id]
  coloc_input[, t2_excluded := t2_uniqueid %in% exclude_dt$unique_id]
  coloc_input <- coloc_input[!(coloc_input$t1_excluded & coloc_input$t2_excluded)]

  message(nrow(coloc_input), " coloc tests remaining after exclusion")
}

if (!is.null(opt$exclude_oneside)) {
  exclude_1side <- readr::read_delim(opt$exclude_oneside, delim = "\t", col_names = TRUE)
  exclude_1side <- as.data.table(exclude_1side)
  
  message("Found ", nrow(exclude_1side), " entries in one side exclusion")
  exclude_1side[, unique_id := paste(study_id, phenotype_id, sep = "_")]

  # Exclude studies where either t1 or t2 are in the exclusion list
  coloc_input[, t1_excluded := t1_uniqueid %in% exclude_1side$unique_id]
  coloc_input[, t2_excluded := t2_uniqueid %in% exclude_1side$unique_id]
  coloc_input <- coloc_input[!(coloc_input$t1_excluded | coloc_input$t2_excluded)]

  message(nrow(coloc_input), " coloc tests remaining after exclusion")
}

# Remove accessory columns
coloc_input[, t1_uniqueid := NULL]
coloc_input[, t2_uniqueid := NULL]
coloc_input[, t1_excluded := NULL]
coloc_input[, t2_excluded := NULL]

# Save coloc guide table
data.table::fwrite(coloc_input, file = opt$output_file)
message("Coloc guide table saved to ", opt$output_file)
