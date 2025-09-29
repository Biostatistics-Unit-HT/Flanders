# Load packages
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(coloc))
suppressMessages(library(susieR))
#suppressMessages(library(bigsnpr))
suppressMessages(library(stringi))
suppressMessages(library(stringr))
#suppressMessages(library(reshape2))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(plyr))
#suppressMessages(library(Gviz))
suppressMessages(library(Matrix))
suppressMessages(library(Rfast))
suppressMessages(library(dplyr))




#### run_dentist #### 
# Preparation of files necessary to perform DENTIST - same files will be used for COJO!!
run_dentist <- function(D=dataset_aligned
                        , locus_chr=opt$chr
                        , locus_start=opt$start
                        , locus_end=opt$end
                        , bfile="ld_panel"
                        , maf.thresh=1e-4
                        , random.number="ZUlGe4EnYqGkubYrApHu"
                        , dentist.bin="DENTIST"
){
  
  
  # Save list of snps included in the locus    
  locus_only.snp <- D |> 
    dplyr::filter(CHR==locus_chr, BP >= locus_start, BP <= locus_end) |>
    dplyr::pull(SNP)
  write(locus_only.snp, ncol=1,file=paste0(random.number,"_locus_only.snp.list"))
  
  # Prepare subset of plink LD files    
  exit_status = system(paste0("plink2 --pfile ", bfile," --extract ",random.number,"_locus_only.snp.list --maf ", maf.thresh, " --make-bed --out ", random.number))
  
  # Raise an error if the external command fails
  if (exit_status != 0) {
    cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
    quit(status = 1, save = "no")
  }

  #### Check if any SNP was left after extracting and filtering!
  snsp_extracted <- system(paste0("grep -E '[0-9]+ variants? remaining after main filters\\.' ", random.number, ".log"), intern = TRUE)
  snsp_extracted <- as.numeric(gsub("(\\d+) variants? remaining after main filters.", "\\1", snsp_extracted))
  
  
  if(length(snsp_extracted) > 0){
    
    # Format gwas sum stat input
    D <- D |>
      dplyr::select("SNP","A1","A2","freq","b","se","p","N","snp_original","type", any_of(c("sdY","s")))
    
    fwrite(D |> dplyr::select(-snp_original, -type, -any_of(c("sdY", "s"))), # to match with input required by Dentist
           file=paste0(random.number,"_sum.txt"), row.names=F,quote=F,sep="\t", na=NA)
    
    # We don't catch the exit status of the system call here, as we want to continue even if DENTIST fails
    system(paste0(dentist.bin, "/DENTIST_1.3.0.0 --gwas-summary ", random.number,"_sum.txt --bfile ", random.number, " --chrID ", locus_chr,  " --extract ", random.number, "_locus_only.snp.list --out ", random.number, " --thread-num 1"))
    
    if (file.exists(paste0(random.number, ".DENTIST.short.txt"))){ ### check that output was produced
      # Remove SNPs pointed out by DENTIST and proceed with COJO
      dentist_exclude <- fread(paste0(random.number, ".DENTIST.short.txt"), data.table = F, header = F) 
      if (nrow(dentist_exclude)>0){ ### check that output produced isn't empty
        locus_only.snp <- setdiff(locus_only.snp, dentist_exclude[,1])
        # SAVE UPDATED LIST OF SNPS!!!        
        write(locus_only.snp, ncol=1,file=paste0(random.number,"_locus_only.snp.list"))
      }
    }
  } else {
    cat(paste0(snsp_extracted, " variants remaining in the LD reference panel after SNPs extraction and MAF filter"))
    system(paste0("rm ", random.number, "*"))
    quit(save = "no", status = 0, runLast = FALSE)  # Exit the script gracefully
  }
}



#### prep_susie_ld ####
# Prepare LD matrix for SUSIE
prep_susie_ld <- function(
    D=dataset_aligned,
    locus_chr=opt$chr,
    locus_start=opt$start,
    locus_end=opt$end,
    bfile=opt$bfile,
    maf.thresh=opt$maf,
    random.number="ZUlGe4EnYqGkubYrApHu",
    skip_dentist=opt$skip_dentist
){
  
  if (skip_dentist == TRUE){
    # Save list of snps included in the locus    
    locus_only.snp <- D |> 
      dplyr::filter(CHR==locus_chr, BP >= locus_start, BP <= locus_end) |>
      dplyr::pull(SNP)
    write(locus_only.snp, ncol=1,file=paste0(random.number,"_locus_only.snp.list"))
  }
  
  ### --export A include-alt --> creates a new fileset, after sample/variant filters have been applied - A: sample-major additive (0/1/2) coding, suitable for loading from R 
  exit_status = system(paste0("plink2 --pfile ", bfile, " --extract ", random.number, "_locus_only.snp.list --maf ", maf.thresh, " --export A include-alt --out ", random.number))
  
  # Check if the command failed
  if (exit_status != 0) {
    # Check if the error is due to no variants remaining
    plink_log <- paste0(random.number, ".log")
    if (file.exists(plink_log)) {
      log_content <- readLines(plink_log)
      if (any(grepl("No variants remaining after main filters", log_content))) {
        cat("Warning: No variants remaining after filtering for locus. Skipping this locus.\n")
        system(paste0("rm ", random.number, "*"))  # Clean up temporary files
        return(NULL)  # Skip further processing for this locus
      }
    }
    # If the error is not due to no variants, stop the pipeline
    stop("Error: External command failed with exit code: ", exit_status)
  }

  geno <- fread(paste0(random.number, ".raw"))[,-c(1:6)] ### First 6 columns are FID, IID, PAT, MAT, SEX and PHENOTYPE
  
  # Check which SNPs have the same genotype for all samples and remove them
  not_same_geno <- which(sapply(geno, function(x) length(unique(x)) > 1))
  geno <- geno[, ..not_same_geno]
  
  # split the SNP names into rsID, effective and other alleles
  snp_info <- strsplit(colnames(geno), "_|\\(/|\\)")
  snp_info <- Reduce(rbind,snp_info) |>
    base::data.frame()

# If there's only one SNP in the plink .raw file, it will need further formatting
  if(ncol(snp_info)==1){
    snp_info <- data.frame(t(snp_info))
  }
  
  colnames(snp_info) <- c('SNP','ea','oa')
  rownames(snp_info) <- NULL
  colnames(geno) <- snp_info$SNP
  
  ##### Ideally to have in the bfile processing step  
  # Remove SNPs with duplicated ids (all occurrencies!)
  dup_snps_index <- which((duplicated(snp_info$SNP) | duplicated(snp_info$SNP, fromLast = TRUE)))
  
  if(length(dup_snps_index)>0){
    snp_info <- snp_info[-dup_snps_index, ]
    geno <- geno[, -..dup_snps_index]
  }
  #####
  
  # check for which columns genotypes should be reverted
  index_to_revert <- which(!(snp_info$ea <= snp_info$oa))
  
  # Switch geno --> from 0 to 2 and vice-versa
  # Function to switch 0 to 2 and 2 to 0
  switch_0_2 <- function(x) {
    switched <- (x*-1)+2
    return(switched)
  }
  
  # Apply the transformation only to specified columns - if index_to_revert is empty, skip this step
  if(length(index_to_revert)>0){
    geno[, (index_to_revert) := lapply(.SD, switch_0_2), .SDcols = index_to_revert]
  }
  
  # Impute missing genotypes with mean value  
  geno <- apply(geno, 2, function(x) {x[is.na(x)] <- mean(x,na.rm=TRUE); return(x)})
  # Correlation matrix
  #ld <- cor(geno) #### NB: don't square it!!!!
  X_scaled <- scale(geno)  # Standardize columns
  ld <- crossprod(X_scaled) / (nrow(geno) - 1) # Same as cor(), but faster
  
  system(paste0("rm ", random.number, "*"))
  return(ld)
}



### hcolo.cojo.ht ###
hcolo.cojo.ht=function(df1 = conditional.dataset1,
                       df2 = conditional.dataset2,
                       p1=1e-4,
                       p2=1e-4,
                       p12=1e-5
                       ){
  
  df1 <- df1 |> dplyr::rename("lABF.df1"="lABF")
  df2 <- df2 |> dplyr::rename("lABF.df2"="lABF")
  
  p1 <- coloc:::adjust_prior(p1, nrow(df1), "1")
  p2 <- coloc:::adjust_prior(p2, nrow(df2), "2")
      
  merged.df <- merge(df1, df2, by = "snp")
  p12 <- coloc:::adjust_prior(p12, nrow(merged.df), "12")
    
  if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")
      
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- coloc:::logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)
      
  pp.abf <- coloc:::combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
      
  colo.res <- list(summary=results, results=merged.df, priors=c(p1=p1,p2=p2,p12=p12))
  class(colo.res) <- c("coloc_abf", class(colo.res))
      
## Save coloc summary        
  colo.sum <- data.frame(t(colo.res$summary))
      
## Save coloc result by SNP
  colo.full_res <- colo.res$results |> dplyr::select(snp,lABF.df1,lABF.df2,SNP.PP.H4)

## Organise all in a list ( composed of summary + results)
  coloc.final <- list(summary=colo.sum, results=colo.full_res)
  return(coloc.final)
}


# Sodbo Sharapov, Human Technopole (c) 2024

#' Filters SuSiE output
#'
#' @description
#' This function filters the output from SuSiE and outputs Susie class object.
#' It return the object of class "susie" with filtered:
#' susie_output$sets$cs
#' susie_output$sets$purity
#' susie_output$sets$cs_index
#' susie_output$sets$coverage
#'
#' @param susie_output A SuSiE_R output (e.g. susieR::fitted_rss()).
#' @param pval A numeric vector representing P-values of the SNPs in the locus, tested in the susie_output.
#' @param cs_bf_thr A numeric value representing the credible set BF threshold for filtering credible sets (default is 3).
#' @param signal_pval_threshold A numeric value representing top SNP P-value threshold for filtering credible sets (default is 1).
#' @param purity_mean_r2_threshold A numeric value representing the threshold for purity mean r2 QC metrics for filtering credible sets (default is 0).
#' @param purity_min_r2_threshold A numeric value representing the threshold for purity min r2 QC metrics for filtering credible sets (default is 0).
#' @param verbose A boolean value. If TRUE, the log information will be printed.
#'
#' @return A Susie object containing filtered fine-mapped credible sets.
#'
#' @details
#' This function takes the output from SuSiE, filters credible sets by purity metrics and CS lbf,
#' and returns a Susie object.
#'
#' @examples
#' \dontrun{
#' susie_filtered <- susie.cs.ht(
#'     susie_output = susie_output,
#'     pval = pval,
#'     cs_bf_thr = 3, 
#'     signal_pval_threshold = 1,
#'     purity_mean_r2_threshold = 0,
#'     purity_min_r2_threshold = 0,
#'     verbose = TRUE
#' )
#' }
#' @export
susie.cs.ht <- function(
    susie_output = NULL,
    pval = NULL,
    cs_bf_thr = 3,
    signal_pval_threshold = 1,
    purity_mean_r2_threshold = 0,
    purity_min_r2_threshold = 0,
    verbose = FALSE
) {
  
  # Head of summary.susie()$cs:
  # cs cs_log10bf cs_avg_r2 cs_min_r2 variable
  # 2  22.803334 1.0000000 1.0000000 371
  # 3  13.595152 1.0000000 1.0000000 255
  # 8   4.890725 1.0000000 1.0000000 492
  # 2.738234 0.9994622 0.9994622 200,202
  
  cs_info <- susieR::summary.susie(susie_output)$cs

  # Get names of the CS that have bf above threshold
  good_bf_cs <- cs_info$cs[cs_info$cs_log10bf >= log10(cs_bf_thr)]

  # Get names of the CS that have mean purity above threshold
  good_purity_mean_cs <-  cs_info$cs[cs_info$cs_avg_r2 >= purity_mean_r2_threshold]
  
  # Get names of the CS that have min purity above threshold
  good_purity_min_cs <-  cs_info$cs[cs_info$cs_min_r2 >= purity_min_r2_threshold]
  
  if(!is.null(pval)) {
    # Get vector of indexes of SNPs per each CS
    min_pval_per_cs <- lapply(
      strsplit(cs_info$variable,split = ","),
      as.numeric)
    
    # Get vector of P-values of top SNP per each CS
    min_pval_per_cs <- sapply(min_pval_per_cs,function(x) min(pval[x]))
    
    # Get names of the CS that have index SNP with P-value below threshold
    good_pval_cs <- cs_info$cs[min_pval_per_cs <= signal_pval_threshold]
    
  } else {
    
    good_pval_cs <- cs_info$cs
    
  }
  
  # Get intersect of CS that satisfy all the thresholds
  cs_to_keep <- Reduce(intersect,list(good_bf_cs, good_purity_mean_cs, good_purity_min_cs,good_pval_cs))
  
  susie_output_filtered <- susie_output
  
  # Print verbose information
  if(verbose){
    message(paste0(c("CS with bf above threshold: ",good_bf_cs),collapse = " "))
    message(paste0(c("CS with mean purity above threshold: ",good_purity_mean_cs),collapse = " "))
    message(paste0(c("CS with min purity threshold: ",good_purity_min_cs),collapse = " "))
    message(paste0(c("P-values of the index SNP for each CS: ",min_pval_per_cs),collapse = " "))
    message(paste0(c("CS with index SNP P-value below threshold: ",good_pval_cs),collapse = " "))
    message(paste0(c("CS which sutisfy all thresholds: ",cs_to_keep),collapse = " "))
  }
  
  # If there is at least one CS that passed threshold - report it, otherwise return NULL
  if(length(cs_to_keep)!=0){
    
    susie_output_filtered$sets$cs <- susie_output_filtered$sets$cs[paste0("L",cs_to_keep)]
    susie_output_filtered$sets$purity <- susie_output_filtered$sets$purity[paste0("L",cs_to_keep),]
    susie_output_filtered$sets$coverage <- susie_output_filtered$sets$coverage[which(susie_output_filtered$sets$cs_index %in% cs_to_keep)]
    susie_output_filtered$sets$cs_index <- cs_to_keep
    
    return(susie_output_filtered)
    
  } else {
    message("None of the credible sets satisfy selected thresholds")
    return(NULL)
  } 
  
}


### From lABF to PP - sent by Nicola
logbf_to_pp_ht = function(bf=all.coloc.join$ABF.tot) {
  denom=coloc:::logsum(bf)
  exp(bf  - denom)
}


### Obtained conditional beta and se from Susie output
get_beta_se_susie <- function(sus,L_index){
  
  se = sqrt(sus$mu2[L_index,] - (sus$mu[L_index,])^2)/sus$X_column_scale_factors
  beta = (sus$mu[L_index,]/sus$X_column_scale_factors)
  
  return(list(beta = beta, se = se))
  
}


# Helper function to run susie_rss with error handling
run_susie_w_tryCatch <- function(
    D_sub,
    D_var_y,
    susie_ld,
    L = L,
    coverage = coverage_value,
    max_iter = max_iter,
    min_abs_corr = NULL
) {  
  tryCatch({
    susie_rss(
      bhat = D_sub$b, 
      shat = D_sub$se, 
      n = max(D_sub$N), 
      R = susie_ld, 
      var_y = D_var_y,
      estimate_residual_variance = FALSE,
      L = L,
      coverage = coverage,
      max_iter = max_iter,
      min_abs_corr = min_abs_corr
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    message(msg)
    if (grepl("The estimated prior variance is unreasonably large", msg)) {
      return("SKIP_TO_L1")
    }
    stop(e)  # Re-throw other errors
  })
}


# Run SuSiE with retries
run_susie_w_retries <- function(
    D_sub,
    D_var_y,
    susie_ld,
    L = L,
    coverage = coverage,
    max_iter = max_iter,
    min_coverage = min_coverage,
    min_abs_corr = NULL
){
  
  set.seed(1) # To ensure reproducibility
  fitted_rss <- NULL  # Initialize
  coverage_value_updated <- coverage
  
  # First run!  
  fitted_rss <- run_susie_w_tryCatch(
    D_sub,
    D_var_y,
    susie_ld,
    L = L,
    coverage = coverage,
    max_iter = max_iter,
    min_abs_corr = 0.5 ## default value
  )
  
  # "Estimated prior variance is unreasonably large" error - jump to L=1
  if (identical(fitted_rss, "SKIP_TO_L1")) {
    message("Re-running fine-mapping with L=1")
    
    fitted_rss <- run_susie_w_tryCatch(
      D_sub,
      D_var_y,
      susie_ld,
      L = 1,
      coverage = coverage,
      max_iter = max_iter,
      min_abs_corr = 0 ## do not filter for anything
    )
    fitted_rss$comment_section <- "The estimated prior variance is unreasonably large. This is usually caused by mismatch between the summary statistics and the LD matrix. Please check the input."
    
    # If cs are still NULL, return NULL
    if(is.null(fitted_rss$sets$cs)){
      return(NULL)  # return early - stop here function!
    }
#    write.table(msg, "failed_susie.txt", row.names = FALSE, col.names = FALSE)
#    quit(save = "no", status = 0, runLast = FALSE)  # Exit the script gracefully # would not work in a loop setting
#    return(NULL)
  } else {
    fitted_rss$comment_section <- NA
  }


  # If result or cs is NULL (but not for L=1 cases) - retry lowering coverage, if still NULL jump to L=1
  while ( (is.null(fitted_rss) || is.null(fitted_rss$sets$cs)) && coverage_value_updated >= min_coverage) {
    coverage_value_updated <- coverage_value_updated - 0.01
    message("Re-running with reduced coverage: ", coverage_value_updated)
    
    fitted_rss <- run_susie_w_tryCatch(
      D_sub,
      D_var_y,
      susie_ld,
      L = L,
      coverage = coverage_value_updated,
      max_iter = max_iter,
      min_abs_corr = 0.5 ## default value
    )
    fitted_rss$comment_section <- NA
  }

  # If still no credible sets, try L = 1 as last resort
  if (coverage_value_updated < min_coverage) {
    fitted_rss <- run_susie_w_tryCatch(
      D_sub,
      D_var_y,
      susie_ld,
      L = 1,
      coverage = coverage,
      max_iter = max_iter,
      min_abs_corr = 0 ## do not filter for anything
    )
    fitted_rss$comment_section <- paste0("Final attempt of fine-mapping failed, reached minimum coverage of ", min_coverage, ". Re-run with L=1")

    # If cs are still NULL, return NULL
    if(is.null(fitted_rss$sets$cs)){
      return(NULL)  # return early - stop here function!
    }
  }
    
    
  # Check if susie converged in given number of iterations - if not, jump to L=1
  if (!(fitted_rss$converged)) {
    message("IBSS algorithm did not converge in ", max_iter, " iterations")
    message("Re-running fine-mapping with L=1")
    fitted_rss <- run_susie_w_tryCatch(
      D_sub,
      D_var_y,
      susie_ld,
      L = 1,
      coverage = coverage,
      max_iter = max_iter,
      min_abs_corr = 0 ## do not filter for anything
    )
    fitted_rss$comment_section <- paste0("IBSS algorithm did not converge in ", max_iter, " iterations! Please check consistency between summary statistics and LD matrix. See https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html")
    
    # If cs are still NULL, return NULL
    if(!(fitted_rss$converged)){
      return(NULL)  # return early - stop here function!
    }
  }
    
  return(fitted_rss)
}


#' Expand SuSiE Credible Sets by Thresholding
#'
#' Expands each credible set (CS) from a fitted \code{susie} or \code{susie_rss} 
#' object by applying a thresholding rule to log Bayes factors (\code{lbf_variable}). 
#' The expanded set is defined as the union of the original CS and 
#' SNPs identified via \code{\link{find_threshold}}.
#'
#' @param fitted A \code{susie} or \code{susie_rss} object, typically the 
#'   output of \code{susieR::susie_rss()}.
#'
#' @return A named list of integer vectors, one per credible set.  
#'   Each vector contains SNP indices for the union of the original and 
#'   expanded credible sets, with SNP IDs as names.
#'
#' @details 
#' The function retrieves original CSs, then applies \code{find_threshold()}
#' to the corresponding lBF values. The expanded CS is the union of both sets.
#'
#' @seealso \code{\link{find_threshold}}
#' @export
#'
#' @examples
#' \dontrun{
#' fitted <- susieR::susie_rss(z, R, n = N)
#' expanded <- expand_cs(fitted)
#' names(expanded)
#' }
expand_cs <- function(fitted) {
    cs_list  <- fitted$sets$cs
    snp_names <- colnames(fitted$lbf_variable)
    
    expanded <- lapply(names(cs_list), function(i) {
        idx_original <- cs_list[[i]]
        names(idx_original) <- snp_names[idx_original]  # name with SNP IDs
        i_num <- as.integer(gsub("L", "", i))

        lbf_vec <- fitted$lbf_variable[i_num, ]
        idx_expanded <- find_threshold(lbf_vec, original_idx = idx_original)
    
        # Union of both
        idx_union <- union(idx_original, idx_expanded)
        names(idx_union) <- snp_names[idx_union]
        
        idx_union
    })
    names(expanded) <- names(cs_list)
    return(expanded)
}



#' Identify Expanded Credible Set SNPs via Thresholding
#'
#' Given a vector of SNP scores (e.g. lABFs), this function determines an 
#' optimal threshold that best separates the values into two groups. 
#' SNPs above the threshold are considered part of the expanded credible set.
#'
#' @param vector Numeric vector of SNP-level scores.
#' @param original_idx Optional integer vector of SNP indices corresponding 
#'   to the original credible set. If no expansion is possible 
#'   (e.g. variance of \code{vector} is zero), these indices are returned instead.
#'
#' @return An integer vector of SNP indices corresponding to the expanded credible set.  
#'   If no expansion is possible, the original indices are returned.
#'
#' @seealso \code{\link{expand_cs}}
#' @export
#'
#' @examples
#' lbf_vec <- rnorm(100)
#' cs_idx <- which(lbf_vec > 2)
#' find_threshold(lbf_vec, original_idx = cs_idx)
find_threshold <- function(vector, original_idx = NULL) {
  if (var(vector) != 0) {
    thres <- optim(
      fn = function(th) thresholder(th, vector),
      p = 0,
      method = "Brent",
      lower = min(vector),
      upper = max(vector)
    )
    idx <- which(vector > thres$par)
    return(idx)
  } else {
    # Return original indices if no expansion is possible
    return(original_idx)
  }
}


#' Objective Function for Threshold Optimization
#'
#' Internal helper used by \code{\link{find_threshold}}. 
#' Given a numeric threshold and a vector of SNP scores, it splits the vector 
#' into two groups and returns the variance of residuals from a simple linear model.
#'
#' @param threshold Numeric. Threshold value used to split the vector into groups.
#' @param vector Numeric vector of SNP scores (e.g. lABFs).
#'
#' @return A single numeric value: the variance of residuals from the fitted model.
#'
#' @note This function is not exported. It is used internally by \code{find_threshold()}.
#'
#' @keywords internal
thresholder <- function(threshold, vector) {
  var(residuals(lm(vector ~ ifelse(vector > threshold, "A", "B"))))
}


#' Transfer SuSiE RSS fine-mapping results into an AnnData object
#'
#' This function takes a list of fine-mapping results from susie and associated metadata,
#' attaches conditional effect size estimates, and organizes the results into an
#' [anndata::AnnData] object with sparse matrices for lABFs, betas, and standard
#' errors. Metadata about credible sets and variants is stored in the `obs` and
#' `var` slots of the output object.
#'
#' @param finemap_list A named list of susie fine-mapping results. Each element is
#'   expected to contain at least:
#'   - `lbf_variable`: a matrix of log Bayes factors
#'   - `sets`: credible set information (`cs_index`, `cs`, `purity`,
#'     `requested_coverage`)
#'   - `metadata`: a data frame with study, phenotype, and variant annotations
#'   - `comment_section`: optional notes from the fine-mapping run.
#' @param cs_indices A named list parallel to `finemap_list`, containing indices
#'   of expanded credible set SNPs to be used instead of the original indices.
#' @param beta_se_cond A named list parallel to `finemap_list`, where each element
#'   contains conditional effect sizes and their standard errors. Each entry should
#'   be a list of lists with elements `beta` and `se`.
#'
#' @details
#' The function performs the following steps for each fine-mapping result:
#' - Identifies and flags credible SNPs according to credible set expansion
#' - Computes and combines metadata about credible sets
#' - Builds sparse matrices of lABFs, betas, and SEs
#' - Constructs an [anndata::AnnData] object with:
#'   - `X`: sparse lABF matrix
#'   - `layers`: beta and SE sparse matrices
#'   - `obs`: credible set metadata
#'   - `var`: variant metadata
#'
#' @return An [anndata::AnnData] object containing fine-mapping summary
#'   statistics and metadata, ready for downstream analysis or visualization.
#'
#' @importFrom dplyr mutate select filter pull
#' @importFrom tidyr unnest_longer
#' @importFrom Matrix sparseMatrix
#' @importFrom anndata AnnData
#' @examples
#' \dontrun{
#' ad <- from_susie_to_anndata(finemap_list, cs_indices, beta_se_cond)
#' }
#'
from_susie_to_anndata <- function(finemap_list=NULL, cs_indices=NULL, beta_se_cond=NULL) {

  lABFs_list <- list()
  min_res_labf_vec <- c()
  top_pvalue_vec <- c()
  purity_df <- data.frame()
  comment_section <- c()
  metadata_df <- data.frame()

  for (finemap_name in names(finemap_list)) {
    finemap <- finemap_list[[finemap_name]]
    cs_index <- cs_indices[[finemap_name]]
    beta_se <- beta_se_cond[[finemap_name]]

    # Compute lABFs
    lABFs <- lapply(finemap$sets$cs_index, function(cs) {
      data.frame(
        SNP = colnames(finemap$lbf_variable),
        lABF = finemap$lbf_variable[cs, ],
        is_cs = FALSE,
        bC = beta_se[[cs]]$beta, # Add conditional beta and se
        bC_se = beta_se[[cs]]$se
      )
    })

    # Get credible SNPs for each CS
    credible_snps <- lapply(finemap$sets$cs_index, function(cs) {
      index <- paste0("L", cs)
      #susie_get_cs(finemap, coverage = finemap$sets$requested_coverage)$cs[[index]]
      cs_index[[index]] ### take indeces of expanded cs rather than original one!
    })

    # Mark credible SNPs
    for (i in seq_along(lABFs)) {
      lABFs[[i]]$is_cs[credible_snps[[i]]] <- TRUE
    }

    # Add info to metadata, collapse study and pheno ID in a single TRAITID
    finemap$metadata <- finemap$metadata |>
      dplyr::mutate(
        snp = list(sapply(lABFs, function(x) x$SNP[which.max(x$lABF)])),
        L_index = list(names(finemap$sets$cs)),
        trait = ifelse(unique(TYPE=="gwas"), study_id, paste0(study_id, ":", phenotype_id))
      ) |>
      tidyr::unnest_longer(c(snp, L_index)) %>%
      dplyr::select(-study_id, -phenotype_id)

    min_res_labf <- sapply(lABFs, function(x) unique(min(x$lABF)))
    top_pvalue <- sapply(lABFs, function(x) min(2 * pnorm(abs(x$bC) / x$bC_se, lower.tail = FALSE)))
    logsum_lABF <- sapply(lABFs, function(x) coloc:::logsum(x$lABF))

    # Store results
#    names(lABFs) <- credible_set_names
    lABFs_list <- c(lABFs_list, lABFs)
    min_res_labf_vec <- c(min_res_labf_vec, min_res_labf)
    top_pvalue_vec <- c(top_pvalue_vec, top_pvalue)
    purity_df <- rbind(purity_df, finemap$sets$purity |> dplyr::mutate(logsum_lABF=logsum_lABF)) 
    comment_section <- c(comment_section, rep(finemap$comment_section, length(finemap$sets$cs_index)))
    comment_section[is.na(comment_section)] <- "NaN"
    metadata_df <- rbind(metadata_df, finemap$metadata) 
  }

#### Create names for cs ####
  cs_names_vec <- paste(
    paste0("chr", metadata_df$chr),
    metadata_df$trait,
    metadata_df$snp,
    metadata_df$L_index,
    sep = "::"
  )

# Prepare `obs_df` metadata
  obs_df <- metadata_df |> dplyr::select(-L_index) |> dplyr::rename(type=TYPE)
  obs_df$top_pvalue <- top_pvalue_vec
  obs_df$min_res_labf <- min_res_labf_vec
#  obs_df$panel <- NaN
  obs_df$cs_name <- cs_names_vec
  names(lABFs_list) <- obs_df$cs_name ### important for extracton of beta and se
  obs_df$a1 <- gsub(".*:(\\w+):(\\w+)$", "\\1", obs_df$snp)
  obs_df$a0 <- gsub(".*:(\\w+):(\\w+)$", "\\2", obs_df$snp)
  obs_df$freq <- NaN

  obs_df <- obs_df |>
    dplyr::mutate(
      bC = map2_dbl(cs_name, snp, ~ lABFs_list[[.x]] |> filter(SNP == .y) %>% pull(bC)),
      bC_se = map2_dbl(cs_name, snp, ~ lABFs_list[[.x]] |> filter(SNP == .y) %>% pull(bC_se))
    )

  obs_df <- cbind(obs_df, purity_df)
  obs_df$comment_section <- comment_section

  # Prepare sparse matrices
  all_snps <- unique(unlist(lapply(lABFs_list, function(x) x$SNP)))
  element_indices <- setNames(seq_along(all_snps), all_snps)

  create_sparse_matrix <- function(value_name, credible_sets=cs_names_vec) {
    Matrix::sparseMatrix(
      i = unlist(lapply(seq_along(credible_sets), function(i) rep(i, length(lABFs_list[[credible_sets[i]]]$SNP)))),
      j = unlist(lapply(credible_sets, function(cs) element_indices[lABFs_list[[cs]]$SNP])),
      x = unlist(lapply(lABFs_list, function(x) x[[value_name]])),
      dims = c(length(credible_sets), length(all_snps)),
      dimnames = list(credible_sets, all_snps)
    )
  }

  lABF_matrix_sparse <- create_sparse_matrix("lABF")
  beta_matrix_sparse <- create_sparse_matrix("bC")
  se_matrix_sparse <- create_sparse_matrix("bC_se")

  # Create AnnData object
  ad <- anndata::AnnData(X = lABF_matrix_sparse)
  ad$layers[["beta"]] <- beta_matrix_sparse
  ad$layers[["se"]] <- se_matrix_sparse
  rownames(obs_df) <- cs_names_vec
  ad$obs <- obs_df

  var_df <- data.frame(
    snp = all_snps,
    chr = gsub(".*(\\d+):(\\d+):\\w+:\\w+", "\\1", all_snps),
    pos = gsub(".*(\\d+):(\\d+):\\w+:\\w+", "\\2", all_snps)
  )
  rownames(var_df) <- var_df$snp
  ad$var <- var_df

  return(ad)
}
