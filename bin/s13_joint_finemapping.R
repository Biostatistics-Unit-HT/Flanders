#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))
suppressMessages(library(Matrix))
suppressMessages(library(flanders))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--annData", default=NULL, help="AnnData cointaing fine-mapping results used for running colocalisation analysis"),
  make_option("--master_table_chunk", default=NULL, help="Chunck of x modules from master table"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read in fine-mapping anndata
ad <- zellkonverter::readH5AD(file = opt$annData, verbose = TRUE, reader = "R")

# Read in master table chunk
master_table_chunk <- fread(opt$master_table_chunk)

# Perform joint fine-mapping
joint_finemapping_obj <- flanders::run_joint_finemapping(ad, master_table_chunk)

message("Joint fine-mapping performed!")

# Remove 0 columns from matrix (any?)
nonzero_cols <- which(colSums(joint_finemapping_obj$joint.lABF.matrix != 0) > 0)
joint_finemapping_obj$joint.lABF.matrix <- joint_finemapping_obj$joint.lABF.matrix[, nonzero_cols, drop = FALSE]

# Save output
# joint_finemapping_obj$master.table  # Access master table with modules, most likely SNP from joint fine-mapping
# joint_finemapping_obj$joint.lABF.matrix # Access sparse matrix collecting joint fine-mapping lABFs (per credible set and SNP)

chunk_n <- gsub(".*_master_table_chunk_(\\d+).tsv", "\\1", opt$master_table_chunk)

fwrite(joint_finemapping_obj$master.table, paste0(opt$coloc_id, "_joint_finemapping_master_table_chunk_", chunk_n, ".tsv"), quote=F, na=NA, sep="\t")
saveRDS(
  list(
    joint.lABF.matrix=joint_finemapping_obj$joint.lABF.matrix,
    cs_metadata=joint_finemapping_obj$cs_metadata,
    master_table_chunk=joint_finemapping_obj$master.table
  ),
  paste0(opt$coloc_id, "_joint_finemapping_lABF_matrix_chunk_", chunk_n, ".rds"))
