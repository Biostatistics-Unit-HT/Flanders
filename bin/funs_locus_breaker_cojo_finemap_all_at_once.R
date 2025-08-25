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
