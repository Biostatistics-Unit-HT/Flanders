#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))
suppressMessages(library(flanders))
suppressMessages(library(dplyr))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--annData", default=NULL, help="AnnData cointaing fine-mapping results used for running colocalisation analysis"),
  make_option("--post_cp_list", default=NULL, help="Path and filenames of R object outputs of run_coloc_percolation()"),
  make_option("--fully_connected_modules", default=NULL, help="Dataframe listing credible sets with fully connected modules and their module label"),
  make_option("--chunk_size", default=NULL, help="Number of modules to analyse together in a chunk"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read in fine-mapping anndata
ad <- zellkonverter::readH5AD(file = opt$annData, verbose = TRUE, reader = "R")

# Read in fully connected modules
fully_connected_modules_df <- fread(opt$fully_connected_modules)

# Read in all run_coloc_percolation() output objects
if(opt$post_cp_list == "empty_post_cp_module_rds.txt"){
  groups <- NULL
  message("Only fully connected modules found - coloc percolation / hyprcoloc_susie not performed")
  
  # Create master table
  master_table <- flanders::create_master_table(
    ad,
    simple_modules = fully_connected_modules_df
  )
  
} else {
  groups_list <- as.vector((read.table(opt$post_cp_list))[,1])
  groups <- lapply(groups_list, readRDS)
  groups <- unlist(groups, recursive = FALSE)
  message("All coloc percolation / hyprcoloc_susie .rds loaded")
  
  # Create master table
  master_table <- flanders::create_master_table(
    ad,
    post_cp_list = groups,
    fully_connected_modules = fully_connected_modules_df
  )
}

message("Master table created")

fwrite(master_table, file = paste0(opt$coloc_id, "_master_table.tsv"), quote=F, na=NA, sep="\t")

# Split master table into n modules
module_groups <- unique(master_table$module)
n_per_file <- as.integer(opt$chunk_size)

# Split module list into chunks of 100
module_chunks <- split(module_groups, ceiling(seq_along(module_groups) / n_per_file))

# Save each chunk to a separate file
for (i in seq_along(module_chunks)) {
  current_modules <- module_chunks[[i]]
  
  master_table_subset <- master_table |> dplyr::filter(module %in% current_modules)
  
  # Save
  fwrite(master_table_subset, file = paste0(opt$coloc_id, "_master_table_chunk_", i, ".tsv"), quote=F, na=NA, sep="\t")
}
