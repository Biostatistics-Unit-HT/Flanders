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
  exit_status = system(paste0("plink2 --bfile ", bfile," --extract ",random.number,"_locus_only.snp.list --maf ", maf.thresh, " --make-bed --out ", random.number))
  
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
  exit_status = system(paste0("plink2 --bfile ", bfile, " --extract ", random.number, "_locus_only.snp.list --maf ", maf.thresh, " --export A include-alt --out ", random.number))
  
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
  if (var(vector) != 0 && !is.na(var(vector))) {
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




#' Extract and Annotate Fine-Mapping Results from SuSiE RSS Output
#'
#' This function processes the output of a \code{susie_rss} fine-mapping run and 
#' returns detailed annotations for each credible set (CS). For each CS, the function 
#' extracts log approximate Bayes factors (lABFs), estimated effect sizes and standard errors,
#' QC metrics (including a log-sum of lABFs), and study metadata. The results are named
#' according to chromosome, study ID, phenotype ID, top SNP, and CS index.
#'
#' @param fitted A \code{susie} or \code{susie_rss} object, typically the output of 
#'   \code{susieR::susie_rss()}.
#' @param D_sub A data.frame containing SNP-level information, with at least the columns:
#'   \code{SNP}, \code{BP} (base-pair position), \code{freq} (allele frequency), and \code{N} (sample size).
#' @param cs_indices A list of credible set indices (e.g. as returned by \code{expand_cs()}), 
#'   where each element is a named vector of SNP indices corresponding to a CS.
#' @param study_id Character. Identifier for the study.
#' @param phenotype_id Character. Identifier for the phenotype.
#' @param chr Integer or character. Chromosome number or label.
#' @param start Integer. Start position of the fine-mapped region.
#' @param end Integer. End position of the fine-mapped region.
#'
#' @return A named list of results, one element per credible set. Each element is itself a list with:
#' \describe{
#'   \item{finemapping_lABFs}{Data frame of SNP-level results, including SNP ID, position, lABF, 
#'         estimated effect size (beta), standard error (bC_se), and whether the SNP belongs to the CS.}
#'   \item{effect}{Data frame describing the top SNP by lABF (alleles, frequency, N, beta, se).}
#'   \item{qc_metrics}{Data frame of purity and coverage metrics for the CS, 
#'         including \code{logsum_lABF}, the log of the summed Bayes factors across SNPs in the set.}
#'   \item{metadata}{Data frame with study and region-level metadata.}
#' }
#'
#' The names of the list elements are constructed as:
#' \code{chr<chr>::<study_id>::<phenotype_id>::<top_snp>::<CS_label>}.
#'
#' @examples
#' \dontrun{
#' fitted <- susieR::susie_rss(z, R, n = N)
#' expanded_cs <- expand_cs(fitted)
#' results <- extract_susie_results(
#'   fitted = fitted,
#'   D_sub = D_sub,
#'   cs_indices = expanded_cs,
#'   study_id = "StudyX",
#'   phenotype_id = "TraitY",
#'   chr = 1,
#'   start = 1e6,
#'   end = 2e6
#' )
#' }
#'
#' @seealso \code{\link{expand_cs}}, \code{\link[susieR]{susie_rss}}
#' @export
extract_susie_results <- function(
    fitted,
    D_sub,
    cs_indices,
    study_id,
    phenotype_id,
    chr,
    start,
    end) {
  
  freq <- setNames(D_sub$freq, D_sub$SNP)
  N    <- setNames(D_sub$N, D_sub$SNP)
  coverage <- fitted$sets$requested_coverage
  
  results <- lapply(fitted$sets$cs_index, function(cs_index) {
    index <- paste0("L", cs_index)
    x <- cs_indices[[index]]  # unified set
    
    beta_se_list <- get_beta_se_susie(fitted, cs_index)
    
    lABF_df <- data.frame(
      SNP   = colnames(fitted$lbf_variable),
      lABF  = fitted$lbf_variable[cs_index, ],
      bC    = beta_se_list$beta,
      bC_se = beta_se_list$se
    )
    
    susie_reformat <- D_sub |>
      dplyr::mutate(is_cs = SNP %in% names(x)) |>
      dplyr::inner_join(lABF_df, by = "SNP") |>
      dplyr::select(SNP, BP, lABF, bC, bC_se, is_cs) |>
      dplyr::rename(snp = SNP, position = BP) |>
      dplyr::arrange(desc(lABF))
    
    # Top SNP
    top <- susie_reformat[1, ]
    snp_top <- top$snp
    chr_pos_a1_a0 <- strsplit(snp_top, ":")[[1]]
    
    effect <- data.frame(
      snp  = snp_top,
      a1   = chr_pos_a1_a0[3],
      a0   = chr_pos_a1_a0[4],
      freq = freq[snp_top],
      N    = N[snp_top],
      beta = top$bC,
      se   = top$bC_se
    )
    
    # QC metrics with log-sum of lABF
    qc_metrics <- fitted$sets$purity[index, ] |>
      dplyr::mutate(
        coverage = coverage,
        L = length(fitted$KL),
        logsum_lABF = coloc:::logsum(susie_reformat$lABF)
      )
    
    metadata_df <- data.frame(
      study_id     = study_id,
      phenotype_id = phenotype_id,
      chr          = chr,
      start        = start,
      end          = end
    )
    
    list(
      finemapping_lABFs = susie_reformat,
      effect            = effect,
      qc_metrics        = qc_metrics,
      metadata          = metadata_df
    )
  })
# ---- Naming each credible set ----
  names(results) <- paste(
    paste0("chr", chr),
    study_id,
    phenotype_id,
    sapply(results, function(x) x$finemapping_lABFs$snp[1]), # top SNP
    names(cs_indices),  # CS label (L1, L2, etc.)
    sep = "::"
  )
  return(results)
}



