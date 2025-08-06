#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(flanders))
suppressMessages(library(zellkonverter))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scRNAseq))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--annData", default=NULL, help="AnnData file with credible set data"),
  make_option("--coloc_guide_table", default=NULL, help="Path and filename of table listing all coloc pairwise tests to perform"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

######  LOAD DATA   ######

# ad <- anndata::read_h5ad(opt$annData)
ad <- sce <- readH5AD(opt$annData, reader = "R")

coloc_guide_table <- fread(opt$coloc_guide_table)

###### COLOCALISTION ######

coloc_res <- anndata2coloc(
  ad,
  coloc_guide_table
  )

fwrite(coloc_res, paste0(#random.number, "_", 
  opt$coloc_id, "_colocalization.table.all.tsv"), quote=F, sep="\t", na=NA)


