#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(igraph))
suppressMessages(library(flanders))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--coloc_results", default=NULL, help="Colocalisation analysis output table reporting pair-wise posterior probabilities results"),
  make_option("--pph4", default=NULL, help="Posterior probability threshold for definifing H4"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

## Read in coloc results
coloc.results <- fread(opt$coloc_results)
coloc.results$PP.H4.abf <- as.numeric(coloc.results$PP.H4.abf)

## Apply function to identify simple modules and check their connectivity
modules <- identify_simple_modules(coloc.results, pp.H4.threshold=as.numeric(opt$pph4), check_connectivity=TRUE)

### Save fully connected modules
fwrite(modules$fully_connected_modules, "modules_fully_connected.tsv", quote = FALSE, na = NA, sep = "\t") # Fully connected

# Split partly connected modules by module and save each to a separate file
split_partly_connected_modules <- split(modules$partly_connected_modules, modules$partly_connected_modules$module)

for(module in names(split_partly_connected_modules)) {
  filename <- paste0("partly_connected_modules_", module, ".tsv")
  fwrite(split_partly_connected_modules[[module]], file = filename, row.names = FALSE, quote=F, na=NA, sep="\t")
}
