#!/usr/bin/env -S Rscript --vanilla

# rds2anndata.R
# This script concatenates multiple RDS files into a single h5ad object,
# fixes the variable table (var) by extracting SNP, chromosome, and position,
# and writes the final AnnData object to an .h5ad file.
#
# Usage:
#   Rscript rds2anndata.R -i <input> -o <output_file>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(anndata))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))



# Define command-line options for the script
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


#' Convert finemapping RDS files into a single AnnData object
#'
#' This function efficiently reads a collection of finemapping `.rds` files and
#' converts them into a single \code{AnnData} object. Each `.rds` file is expected
#' to contain a list of finemapping results (typically 1–10 per file), where each
#' element represents a credible set with associated metadata.
#'
#' For each credible set, the function:
#' \itemize{
#'   \item Extracts SNP-level statistics from \code{finemapping_lABFs}
#'   \item Filters SNPs belonging to the credible set (\code{is_cs == TRUE})
#'   \item Builds sparse matrices for \code{lABF} (stored in \code{X}),
#'         \code{beta}, and \code{se} (stored as layers)
#'   \item Collects credible-set-level metadata into \code{obs}, including:
#'     \itemize{
#'       \item Genomic coordinates and study/phenotype identifiers
#'       \item Minimum p-value within the credible set
#'       \item Minimum lABF across all SNPs
#'       \item Optional effect-level fields (e.g. \code{snp, a1, a0, freq, N, beta, se})
#'       \item Optional QC metrics (e.g. \code{coverage, L, logsum_lABF})
#'     }
#'   \item Collects SNP-level metadata into \code{var} (SNP, chromosome, position)
#' }
#'
#' The function is optimized for large datasets by:
#' \itemize{
#'   \item Reading files sequentially (avoiding full preloading into memory)
#'   \item Using sparse matrix triplet construction (\code{i, j, x})
#'   \item Leveraging \code{data.table} for fast data manipulation
#' }
#'
#' @param finemap_files Character vector of paths to `.rds` finemapping files.
#'   Each file should contain a list of objects with components:
#'   \code{finemapping_lABFs}, \code{metadata}, and optionally \code{effect}
#'   and \code{qc_metrics}.
#' @param panel Character string specifying the SNP genotyping/imputation panel
#'   (stored in \code{ad$obs$panel}). Default is \code{"HRC"}.
#'
#' @return An \code{AnnData} object with:
#' \itemize{
#'   \item \code{X}: sparse matrix of lABF values (credible sets × SNPs)
#'   \item \code{layers[["beta"]]}: sparse matrix of effect sizes
#'   \item \code{layers[["se"]]}: sparse matrix of standard errors
#'   \item \code{obs}: data.frame of credible-set-level metadata
#'   \item \code{var}: data.frame of SNP-level metadata
#' }
#'
#' @details
#' Each element inside each `.rds` file is treated as an independent credible set.
#' Row names of \code{obs} correspond to credible set names extracted from the
#' list names of the `.rds` object. Column names of \code{X} correspond to the
#' union of all SNPs observed across all files.
#'
#' Missing or malformed files are skipped and reported.
#'
#' @examples
#' \dontrun{
#' files <- readLines("input_files.txt")
#' ad <- finemap2anndata_fast(files, panel = "HRC")
#' anndata::write_h5ad(ad, "output.h5ad")
#' }
#'
#' @import data.table
#' @import Matrix
#' @import anndata
#' @export
#'
finemap2anndata_fast <- function(
    finemap_files,
    panel = "HRC"
) {
  
  # Set up lists used to collect AnnData obs/var metadata and sparse matrix triplets  
  obs_list <- list()
  var_list <- list()
  
  i_list <- list()
  snp_list <- list()
  labf_list <- list()
  beta_list <- list()
  se_list <- list()
  
  # Track credible set names and files that failed to load  
  cs_names <- character()
  failed_files <- character()
  
  row_id <- 0L
  
  message("Reading RDS files...")
  
  # Loop over RDS files one at a time to avoid loading all files into memory (rather than kloading in them all at once)
  for (file_idx in seq_along(finemap_files)) {
    f <- finemap_files[file_idx]
    
    # Read one RDS file. If reading fails, record the file and continue
    rds_obj <- tryCatch(
      readRDS(f),
      error = function(e) {
        failed_files <<- c(failed_files, basename(f))
        NULL
      }
    )
    
    if (is.null(rds_obj)) next
    
    # Each inner object inside each RDS file becomes one row in AnnData  
    for (inner_idx in seq_along(rds_obj)) {
      obj <- rds_obj[[inner_idx]]
      
      # Convert finemapping results and metadata to data.table for fast filtering
      fm <- data.table::as.data.table(obj$finemapping_lABFs)
      meta <- data.table::as.data.table(obj$metadata)

      row_id <- row_id + 1L
      
      # Extract credible-set name from rds object name
      cs_name <- names(rds_obj)[inner_idx]
      cs_names[row_id] <- cs_name
      
      # Keep only SNPs belonging to the credible set
      fm_cs <- fm[is_cs == TRUE]
      
      # Store sparse matrix triplets
      # i = AnnData row index, SNP names are converted to column indices later
      i_list[[row_id]] <- rep.int(row_id, nrow(fm_cs))
      snp_list[[row_id]] <- fm_cs$snp
      labf_list[[row_id]] <- fm_cs$lABF
      beta_list[[row_id]] <- fm_cs$bC
      se_list[[row_id]] <- fm_cs$bC_se
      
      # Compute the smallest conditional p-value among credible-set SNPs
      p_values <- stats::pchisq(
        (fm_cs$bC / fm_cs$bC_se)^2,
        df = 1,
        lower.tail = FALSE
      )
      top_pvalue <- min(p_values, na.rm = TRUE)
      
      # Extract per-credible-set effect metadata
      effect_obs <- if (!is.null(obj$effect)) {
        data.table::as.data.table(obj$effect)
      } else {
        data.table::data.table()
      }
      
      # Extract per-credible-set QC metrics
      qc_obs <- if (!is.null(obj$qc_metrics)) {
        data.table::as.data.table(obj$qc_metrics)
      } else {
        data.table::data.table()
      }
      
      # Build obs row: credible-set-level metadata
      obs_list[[row_id]] <- data.table::data.table(
        chr = paste0("chr", meta$chr[1]),
        start = as.numeric(meta$start[1]),
        end = as.numeric(meta$end[1]),
        study_id = meta$study_id[1],
        phenotype_id = meta$phenotype_id[1],
        top_pvalue = top_pvalue,
        min_res_labf = min(fm$lABF, na.rm = TRUE),
        panel = panel,
        cs_name = cs_name
      )[
        ,
        cbind(.SD, effect_obs, qc_obs)
      ]
      
      # Build var metadata: SNP-level information
      var_list[[row_id]] <- fm[, .(
        snp,
        chr = paste0("chr", meta$chr[1]),
        position
      )]
    }
    
    # Progress message every 100 files
    if (file_idx %% 100 == 0) {
      message("Finished ", file_idx, " of ", length(finemap_files), " files")
    }
  }
  
  message("Combining metadata...")
  
  # Combine collected obs and var metadata
  obs_dt <- data.table::rbindlist(obs_list, fill = TRUE)
  var_dt <- unique(
    data.table::rbindlist(var_list, fill = TRUE),
    by = "snp"
  )
  
  # Create SNP-to-column-index mapping for sparse matrix construction
  all_snps <- var_dt$snp
  snp_index <- stats::setNames(seq_along(all_snps), all_snps)
  
  # Flatten triplet lists into vectors
  i_vec <- unlist(i_list, use.names = FALSE)
  snp_vec <- unlist(snp_list, use.names = FALSE)
  j_vec <- unname(snp_index[snp_vec])
  
  labf_vec <- unlist(labf_list, use.names = FALSE)
  beta_vec <- unlist(beta_list, use.names = FALSE)
  se_vec <- unlist(se_list, use.names = FALSE)
  
  message("Building sparse matrices...")
  
  # Matrix dimensions: rows = credible sets, columns = unique SNPs
  dims <- c(length(cs_names), length(all_snps))
  dimnames <- list(cs_names, all_snps)
  
  # Main AnnData matrix: lABF values
  X <- Matrix::sparseMatrix(
    i = i_vec,
    j = j_vec,
    x = labf_vec,
    dims = dims,
    dimnames = dimnames
  )
  
  # AnnData beta layer
  beta <- Matrix::sparseMatrix(
    i = i_vec,
    j = j_vec,
    x = beta_vec,
    dims = dims,
    dimnames = dimnames
  )
  
  # AnnData se layer
  se <- Matrix::sparseMatrix(
    i = i_vec,
    j = j_vec,
    x = se_vec,
    dims = dims,
    dimnames = dimnames
  )
  
  obs_df <- as.data.frame(obs_dt)
  rownames(obs_df) <- cs_names
  
  var_df <- as.data.frame(var_dt)
  rownames(var_df) <- var_df$snp
  
  message("Creating AnnData...")
  
  # Create AnnData with X, obs and var
  ad <- anndata::AnnData(
    X = X,
    obs = obs_df,
    var = var_df
  )
  
  # Add additional sparse matrices as AnnData layers
  ad$layers[["beta"]] <- beta
  ad$layers[["se"]] <- se
  
  message("Done.")
  message("Credible sets: ", length(cs_names))
  message("SNPs: ", length(all_snps))
  message("Failed files: ", length(failed_files))
  
  ad
}



# Read list of finemapped files
input_files <- readLines(opt$input)

# Collect all info into an AnnData
#start <- Sys.time()
ad <- finemap2anndata_fast(
  finemap_files = input_files,
  panel = "HRC"
)
#end <- Sys.time()

# Save anndata
anndata::write_h5ad(ad, filename = opt$output_file)
