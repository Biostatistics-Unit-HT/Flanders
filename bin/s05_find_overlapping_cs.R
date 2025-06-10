#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--coloc_info_table", default=NULL, help="Path and filename of master coloc table produced by individual traits pre-processing"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis"),
  make_option("--chr_cs", default=NULL, help="Chromosomes to split the analysis in")#,
#  make_option("--chunk_size", default=NULL, help="Row number of each chunk the output file will be broken into")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))


# List all cs variants and filter by chr
tb <- readr::read_delim(opt$coloc_info_table, data.table=F)
tb <- as.data.table(tb) |> dplyr::filter(chr==opt$chr_cs)
nrows_tb <- nrow(tb)
tb <- tb  |> dplyr::mutate(cs_name=paste0("cs", seq(1,nrows_tb)))

cs_list <- strsplit(tb$credible_set_snps, ",", fixed=TRUE) # If TRUE match split exactly, otherwise use regular expressions - hoping is more memory efficient! https://stackoverflow.com/questions/55919893/more-memory-efficient-way-than-strsplit-to-split-a-string-into-two-in-r
names(cs_list) <- tb$cs_name

# Get unique elements present in all text vectors
all_elements <- unique(unlist(cs_list))

# Check that there's at list a case of overlapping! Otherwise no point in performing all this - FOR THE MOMENT. Still useful to add it to the massive matrix? But will decide later

if(length(unlist(cs_list)) > length(all_elements)){
  
  element_indices <- setNames(seq_along(all_elements), all_elements)
  
  # Create a sparse matrix with zeros
  matrix_sparse <- Matrix(
    0,
    nrow = length(cs_list),
    ncol = length(all_elements),
    sparse = TRUE,
    dimnames = list(names(cs_list), all_elements)
  )
  
  
  # Create a matrix of indices for each row
  row_indices <- lapply(cs_list, function(vec) {
    element_indices[vec]
  })
  
  # Use this (alternative) - DOESNT WORK --> values must be length 97068, but FUN(X[[1]]) result is length 37
  #row_indices <- vapply(cs_list, function(vec) { element_indices[vec] },
  #                      FUN.VALUE = integer(length(all_elements))
  #)
  
  matrix_sparse[cbind(rep(seq_along(cs_list), lengths(row_indices)), unlist(row_indices))] <- 1
  
  
  
  # Fill in the matrix with 1 where the text vector contains the element
  
  # Save matrix -- GIVE BETTER NAME
  #writeMM(matrix_sparse, "matrix.mtx")
  #system("gzip -f matrix.mtx")
  #####m <- readMM("matrix.mtx.gz") ### To read the gzipped matrix back into R. So you remember
  
  # To get the number of non-zero entries in each column:
  # diff(matrix_table_sparse@p)
  
  # Multiply the transpose of the matrix with itself
  #start <- Sys.time()
  matrix_product <- matrix_sparse %*% t(matrix_sparse)
  #end <- Sys.time()
  #end-start
  
  # Set the diagonal to zero (as each vector will have 1s on diagonal)
  diag(matrix_product) <- 0
  
  shared_elements <- which(matrix_product > 0, arr.ind = TRUE)
  shared_elements_df <- as.data.frame(shared_elements)
  
  # Map row and column indices to credible set names
  shared_elements_df$row <- rownames(shared_elements)
  shared_elements_df$col <- colnames(matrix_product)[shared_elements_df$col]
  
  # Keep unique pairs
  shared_elements_unique <- unique(t(apply(shared_elements_df, 1, sort)))
  colnames(shared_elements_unique) <- c("t1", "t2")
  
  # Retrieve trait info
  coloc_combo <- merge(
    shared_elements_unique,
    tb |> dplyr::select(credible_set_name, study_id, phenotype_id, top_pvalue, path_rds, cs_name),
    by.x = "t1",
    by.y = "cs_name"
  )
  
  coloc_combo <- merge(
    coloc_combo,
    tb |> dplyr::select(credible_set_name, study_id, phenotype_id, top_pvalue, path_rds, cs_name),
    by.x = "t2",
    by.y = "cs_name",
    suffixes = c("_t1", "_t2")
  )
  
  # Remove pair testing different conditional dataset for the same trait (study_id + phenotype_id)
  coloc_combo <- coloc_combo |>
    dplyr::filter(study_id_t1 != study_id_t2 | (study_id_t1==study_id_t2 & phenotype_id_t1 != phenotype_id_t2)) |>
    dplyr::select(-t1, -t2)
  
  colnames(coloc_combo) <- c(
    "t1_credible_set_name", "t1_study_id", "t1_phenotype_id", "t1_top_pvalue", "t1_path_rds",
    "t2_credible_set_name", "t2_study_id", "t2_phenotype_id", "t2_top_pvalue", "t2_path_rds")
  fwrite(coloc_combo, paste0(opt$coloc_id, "_chr", opt$chr_cs, "_coloc_pairwise_guide_table.tsv"), quote=F, na=NA, sep="\t")
}
