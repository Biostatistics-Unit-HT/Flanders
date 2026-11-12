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
  make_option("--fully_connected_modules", default=NULL, help="Dataframe listing credible sets with fully connected modules and their module label"),
  make_option("--post_cp_dir", default=NULL, help="Directory containing module-refinement RDS results"),
  make_option("--chunk_size", default=NULL, help="Number of modules to analyse together in a chunk"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read in fine-mapping anndata
ad <- zellkonverter::readH5AD(file = opt$annData, verbose = TRUE, reader = "R")

# Read in fully connected modules
fully_connected_modules_df <- fread(opt$fully_connected_modules)

# Load-in module-refinement .rds at the following directory
groups_list <- list.files(
  path = opt$post_cp_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

if (length(groups_list) == 0L) {
  message("No refinement results - creating the master table without refined modules")
  
  # Create master table
  master_table <- flanders::create_master_table(
    ad,
    simple_modules = fully_connected_modules_df
  )
  
} else {
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
