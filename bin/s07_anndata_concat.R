#!/usr/bin/env Rscript
# anndata_concat.R
# This script concatenates multiple AnnData files into a single object,
# fixes the variable table (var) by extracting SNP, chromosome, and position,
# and writes the final AnnData object to an .h5ad file.
#
# Usage:
#   Rscript anndata_concat.R -i <input> -o <output_file>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(anndata))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))



# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "List of .h5ad files",
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


#' Convert finemap files to AnnData
#'
#' This function reads finemap `.rds` files, filters rows where `is_cs` is `TRUE`,
#' and creates a sparse matrix of lABF values for SNPs across credible sets.
#' The function then converts the sparse matrix to an AnnData object and
#' includes metadata about the start and end positions of each credible set.
#'
#' @param finemap_files A character vector of file paths to the finemap `.rds` files.
#' @param study_id A character vector specifying the study_id for each of the finemap_files
#' @param snp_panel A character vector specifying the list of SNPs for anndata$X vars. For example,
#' snp_panel can have length of 9M SNPs and resulted AnnData will have X with n(vars) = 9M
#' @param phenotype_id A character vector specifying the phenotype_id for each of the finemap_files
#' @param panel A character string specifying the SNP genotyping/imputation panel
#'
#'
#' @return The AnnData object containing the sparse matrix of lABF values and metadata.
#' @export
#' @import data.table
#' @import dplyr
#' @import anndata
#'
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(anndata)
#' library(data.table)
#' library(dplyr)
#' library(Matrix)
#'

finemap2anndata <- function(
    finemap_files,
    preloaded_list = FALSE,
    study_id = NULL,
    phenotype_id = NULL,
    snp_panel = NULL,
    panel = NULL
){
  
  # Initialize a list to store the filtered data
  filtered_data_list <- list()
  
  # Loop through each file and filter data
  print("Reading finemap files...")
  
  min_res_labf <- c()
  
  # Initialize counters and lists
  failed_files <- c()  # to store failed files
  success_count <- 0      # to count successful reads
  failed_count <- 0       # to count failed reads
  
  # all_snps <- c() # to collect all SNPs, tested in finemapping
  
  # Extract start and end positions from the credible set names
  get_chr_start_end <- function(file_name) {
    # Modify the pattern to capture both 'susie' and 'cojo'
    pattern <- "locus_(chr[0-9XY]+)_([0-9]+)_([0-9]+)(_susie|_cojo)?_finemap\\.rds$"
    
    # Use regmatches and regexec to extract matches
    matches <- regmatches(file_name, regexec(pattern, file_name))
    
    # Extract chromosome, start, and end positions
    chr_start_end <- unlist(lapply(matches, function(x) x[2:4]))
    
    return(chr_start_end)
  }
  
  # List of data.tables with snp, chr and pos columns
  snp_chr_pos <- list()
  
  #setkey(snp_chr_pos,snp)
  
  effect_df <- data.frame()
  qc_metrics_df <- data.frame()
  
  if(!preloaded_list)
    names(finemap_files) <- finemap_files
  
  for (finemap_file in names(finemap_files)) {
    
    # Try to read the .rds file and catch any errors or warnings
    result <- tryCatch({
      if(preloaded_list){
        finemap <- finemap_files[[finemap_file]]
      } else{
        finemap <- readRDS(finemap_file)
      }
      
      # "finemapping_lABFs" "effect"            "qc_metrics"
      
      if(class(finemap) == "list"){
        if(!is.null(finemap$effect)){
          effect_df <- bind_rows(effect_df, finemap$effect)
        } else{
          effect_df <- rbind(effect_df, rep(NA,7))
        }
        if(!is.null(finemap$qc_metrics)){
          qc_metrics_df <- bind_rows(qc_metrics_df, finemap$qc_metrics)
        }else {
          qc_metrics_df <- rbind(qc_metrics_df, rep(NA,3))
        }
        finemap <- finemap$finemapping_lABFs
      }
      finemap <- data.table(finemap)
      if("pC" %in% colnames(finemap) & !("p" %in% colnames(finemap)))
        finemap <- finemap %>% rename(p=pC)
      
      if("bC" %in% colnames(finemap) & !("b" %in% colnames(finemap)))
        finemap <- finemap %>% rename(b=bC)
      
      success_count <- success_count + 1
      
      # If successful, continue with the rest of the operations
      min_res_labf <- c(min_res_labf, min(finemap$lABF))
      
      new_rows <- data.table(
        snp = finemap$snp,
        chr = get_chr_start_end(finemap_file)[1],
        pos = finemap$pos
      )
      
      snp_chr_pos[[length(snp_chr_pos) + 1]] <- new_rows
      
      # Append new rows
      #snp_chr_pos <- rbindlist(list(snp_chr_pos, new_rows), use.names = TRUE)
      
      #snp_chr_pos <- bind_rows(
      #  snp_chr_pos,
      #  data.table(
      #    snp = finemap$snp,
      #    chr = get_chr_start_end(finemap_file)[1],
      #    pos = finemap$pos
      #  )
      #)
      
      #snp_chr_pos <- snp_chr_pos %>% filter(!duplicated(snp))
      
      # Filter rows where is_cs is TRUE
      filtered_finemap <- finemap[is_cs == TRUE]
      
      # Append the filtered data to the list
      filtered_data_list[[basename(finemap_file)]] <- filtered_finemap
      
      NULL  # return NULL on success so tryCatch does not return a value
    }, error = function(e) {
      # On error, append the file name to the failed files list
      print(e)
      failed_files <<- c(failed_files,basename(finemap_file))
      failed_count <<- failed_count + 1
      #message(paste("Error in reading file:", basename(finemap_file), "- Skipping"))
      NULL  # return NULL on error
    })
    n = which(names(finemap_files) == finemap_file)
    if(n %% 100 == 0){
      cat("\rFinished", n, "of", length(finemap_files))
    }
  }
  
  snp_chr_pos <- rbindlist(snp_chr_pos, use.names = TRUE)
  
  snp_chr_pos <- snp_chr_pos %>% filter(!duplicated(snp))
  
  print("")
  print(paste0("We have ",nrow(snp_chr_pos)," SNPs in ", length(filtered_data_list), " credible sets."))
  
  if(preloaded_list)
    finemap_files = names(finemap_files)
  
  if(!is.null(failed_files)){
    study_id <- study_id[-sapply(failed_files,function(x) grep(x,finemap_files))]
    phenotype_id <- phenotype_id[-sapply(failed_files,function(x) grep(x,finemap_files))]
  }
  
  # Output how many files were successfully read and how many failed
  cat("\nNumber of successfully read files:", success_count, "\n")
  cat("Number of failed files:", failed_count, "\n")
  if (!is.null(failed_files)) {
    cat("Failed files:", paste(unlist(failed_files), collapse = ", "), "\n")
  } else {
    cat("Failed files: None\n")
  }
  
  # Collect all unique SNPs
  all_snps <- snp_chr_pos$snp
  
  element_indices <- stats::setNames(seq_along(all_snps), all_snps)
  
  # Collect all credible set names
  credible_sets <- names(filtered_data_list)
  
  # Create a sparse matrix to store lABF values
  lABF_matrix_sparse <- Matrix::Matrix(
    0,
    nrow = length(credible_sets),
    ncol = length(all_snps),
    sparse = TRUE,
    dimnames = list(credible_sets, all_snps)
  )
  
  # Create a sparse matrix to store betas
  beta_matrix_sparse <- Matrix::Matrix(
    0,
    nrow = length(credible_sets),
    ncol = length(all_snps),
    sparse = TRUE,
    dimnames = list(credible_sets, all_snps)
  )
  
  # Create a sparse matrix to store se
  se_matrix_sparse <- Matrix::Matrix(
    0,
    nrow = length(credible_sets),
    ncol = length(all_snps),
    sparse = TRUE,
    dimnames = list(credible_sets, all_snps)
  )
  
  print("Populating sparse matrix...")
  
  top_pvalue <- c()
  
  # Fill the sparse matrix with lABF values
  for (credible_set in credible_sets) {
    credible_data <- filtered_data_list[[credible_set]]
    row_index <- match(credible_set, credible_sets)
    col_indices <- element_indices[credible_data$snp]
    
    lABF_values <- credible_data$lABF
    beta_values <- credible_data$b
    
    if(all(c("b","p") %in% colnames(credible_data)))
      se_values <- abs(
        1/(
          sqrt(
            qchisq(credible_data$p,df=1,lower.tail=FALSE)
          ) / credible_data$b
        )
      )
    
    lABF_matrix_sparse[row_index, col_indices] <- lABF_values
    beta_matrix_sparse[row_index, col_indices] <- beta_values
    se_matrix_sparse[row_index, col_indices] <- se_values
    
    top_pvalue <- c(top_pvalue,min(credible_data$p))
    n = which(credible_sets == credible_set)
    if(n %% 10 == 0){
      cat("\rFinished", n, "of", length(credible_sets))
    }
  }
  
  # This needs to be refactored as now we store chr pos in the ad$var. And it
  # should be properly filled for this "null" snp_panel_matrix_sparse which is
  # appended to the lABF_matrix_sparse
  if(!is.null(snp_panel)){
    
    snp_panel_matrix_sparse <- Matrix::Matrix(
      0,
      nrow = length(credible_sets),
      ncol = length(snp_panel[!snp_panel %in% all_snps]),
      sparse = TRUE,
      dimnames = list(credible_sets, snp_panel[!snp_panel %in% all_snps])
    )
    
    lABF_matrix_sparse <- cbind(lABF_matrix_sparse, snp_panel_matrix_sparse)
    beta_matrix_sparse <- cbind(beta_matrix_sparse, snp_panel_matrix_sparse)
    se_matrix_sparse <- cbind(se_matrix_sparse, snp_panel_matrix_sparse)
    
  }
  
  print("Creating AnnData object...")
  
  # Convert the sparse matrix to an AnnData object
  # TODO - AnnData shoould be created in the end using single call of anndata
  # function - so X, obs, var and layer goes together
  
  ad <- anndata::AnnData(X = lABF_matrix_sparse)
  
  ad$layers[["beta"]] <- beta_matrix_sparse
  ad$layers[["se"]] <- se_matrix_sparse
  
  print("Creating obs meta data...")
  
  # Collect chromosome, start, and end positions for each credible set
  chr_start_end_positions <- t(sapply(credible_sets, get_chr_start_end))
  colnames(chr_start_end_positions) <- c("chr", "start", "end")
  
  # Fill the ad$obs matrix which describes the credible sets
  obs_df <- as.data.frame(chr_start_end_positions, stringsAsFactors = FALSE)
  obs_df$study_id <- study_id
  obs_df$phenotype_id <- phenotype_id
  obs_df$start <- as.numeric(obs_df$start)
  obs_df$end <- as.numeric(obs_df$end)
  obs_df$top_pvalue <- top_pvalue
  obs_df$min_res_labf <- min_res_labf
  obs_df$panel <- panel
  obs_df$cs_name <- credible_sets
  if(nrow(effect_df) > 0 & nrow(qc_metrics_df) > 0){
    obs_df <- bind_cols(obs_df,effect_df,qc_metrics_df)
  }else{
    obs_df <- obs_df
  }
  rownames(obs_df) <- credible_sets # THIS IS VERY IMPORTANT TODO
  
  # Assign the data frame to ad$obs
  ad$obs <- obs_df
  
  print("Creating var meta data...")
  
  # Fill the ad$var matrix which describes the SNPs
  var_df <- snp_chr_pos %>% as.data.frame(, stringsAsFactors = FALSE)
  rownames(var_df) <- snp_chr_pos$snp # THIS IS VERY IMPORTANT TODO
  
  # Assign the data frame to ad$var
  ad$var <- var_df
  
  #print("Writing AnnData to a disk...")
  
  # Save the AnnData object to .h5ad file
  #ad$write_h5ad(output_file)
  
  print("Done...")
  
  return(ad)
}



# List all finemap files
input_files <- strsplit(opt$input, ",")[[1]]



# Function to extract study_id and phenotype_id from filename
get_study_phenotype_id <- function(file_name) {
  pattern <- "^(.*)_chr\\d+:\\d+:\\w+:\\w+_locus_chr\\d+_\\d+_\\d+_.*finemap.rds$"
  matches <- regmatches(file_name, regexec(pattern, file_name))
  study_phenotype_id <- unlist(lapply(matches, function(x) x[2:3]))
  return(study_phenotype_id)
}

# Collect study_id and phenotype_id for each credible set
study_phenotype_ids <- t(sapply(basename(input_files), get_study_phenotype_id))
#' colnames(study_phenotype_ids) <- c("study_id", "phenotype_id")
#' ad_chr_22 <- finemap2anndata(
#'   finemap_files = finemap_chr22_files,
#'   study_id = study_phenotype_ids$study_id,
#'   phenotype_id = study_phenotype_ids$phenotype_id,
#'   output_file = output_file_chr22,
#'   panel = "HRC"
#' )
#'
#' finemap_chr21_files <- list.files(finemap_folder, pattern = "chr21.*\\.rds", full.names = TRUE)
#' finemap_chr21_files <- finemap_chr21_files[!grepl("GWAS", finemap_chr21_files)]
#'
#' # List all ATAC files matching the pattern and filter out those containing "GWAS"
#' finemap_chr22_ATAC_files <- list.files(finemap_folder, pattern = "^ATAC.*\\chr22.*\\.rds", full.names = TRUE)
#' finemap_chr22_ATAC_files <- finemap_chr22_ATAC_files[!grepl("GWAS", finemap_chr22_ATAC_files)]
#'
#' # List all ChipSeq files matching the pattern and filter out those containing "GWAS"
#' finemap_chr22_ChipSeq_files <- list.files(finemap_folder, pattern = "^ChipSeq.*\\chr22.*\\.rds", full.names = TRUE)
#' finemap_chr22_ChipSeq_files <- finemap_chr22_ChipSeq_files[!grepl("GWAS", finemap_chr22_ChipSeq_files)]
#'
#' ad_chr_22 <- finemap2anndata(
#'   finemap_files = finemap_chr22_files,
#'   panel = "HRC"
#' )
#'
#' ad_chr_21 <- finemap2anndata(
#'   finemap_files = finemap_chr21_files,
#'   panel = "HRC"
#' )
#'
#' ad_chr_22_ATAC <- finemap2anndata(
#'   finemap_files = finemap_chr22_ATAC_files,
#'   panel = "HRC"
#' )
#'
#' ad_chr_22_ChipSeq <- finemap2anndata(
#'   finemap_files = finemap_chr22_ChipSeq_files,
#'   panel = "HRC"
#' )
#'
#' chr_21_chr_22_outer <- concat(list(ad_chr_21, ad_chr_22), join = "outer")
#' chr_22_ATAC_chr_22_ChipSeq_outer <- concat(list(ad_chr_22_ATAC, ad_chr_22_ChipSeq), join = "outer")
#' }




library(stringr)

chr22_molQTL_ad$obs$study_id <- str_extract(chr22_molQTL_ad$obs$cs_name, "^[A-Za-z]+_chr[0-9]+")

chr22_molQTL_ad$obs$phenotype_id <- str_match(chr22_molQTL_ad$obs$cs_name, "chr[0-9]+_([^_]+_[0-9]+|ENSG[0-9]+)")[,2]

# Perform colocalization analysis on the GWAS AnnData object
gwas_ad.coloc <- anndata2coloc(gwas_ad)

# View summary of colocalization results
print(gwas_ad.coloc$summary)

# View detailed results by SNP
print(gwas_ad.coloc$results)



