suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--chr", default=NULL, help="Locus chromosome"),
  make_option("--start", default=NULL, help="Locus starting position"),
  make_option("--end", default=NULL, help="Locus ending position"),
  make_option("--phenotype_id", default=NULL, help="Trait for which the locus boundaries have been identified - relevant in cases of molQTLs"),
  make_option("--dataset_aligned", default=NULL, help="GENOME-WIDE munged and aligned dataset file"),
  make_option("--maf", default=1e-04, help="MAF filter", metavar="character"),
  make_option("--bfile", default=NULL, help="Path and prefix name of custom LD bfiles (PLINK format .bed .bim .fam)"),
  make_option("--skip_dentist", default=TRUE, help="Whether to skip the match of SNPs LD between GWAS sum stat and LD reference (performed by DENTIST), and consequent removal of mismatched SNPs"),
  make_option("--cs_thresh", default=0.99, help="Percentage of credible set"),
  make_option("--results_path", default=NULL, help="Path to \"/results\" folder"),
  make_option("--study_id", default=NULL, help="Id of the study")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))


locus_name <- paste0(opt$chr, "_", opt$start, "_", opt$end)
opt$chr <- as.numeric(opt$chr)
opt$start <- as.numeric(opt$start)
opt$end <- as.numeric(opt$end)

# GWAS input
dataset_aligned <- fread(cmd=paste0("tabix ", opt$dataset_aligned, " ", opt$phenotype_id))
colnames(dataset_aligned) <- c("phenotype_id", "snp_original","SNP","CHR","BP","A1","A2","freq","b","se","p","N", "type","temp")

if(unique(dataset_aligned$type=="quant")){
  dataset_aligned <- dataset_aligned %>% dplyr::rename(sdY=temp)
} else {
  dataset_aligned <- dataset_aligned %>% dplyr::rename(s=temp)
}


################################
# Conditional analysis with SUSIE
################################

random.number <- stri_rand_strings(n=1, length=20, pattern="[A-Za-z0-9]")

### If required, run DENTIST to identify mismatches between GWAS sum stats and LD panel
if (opt$skip_dentist == TRUE){
  cat(paste0("I assume you provided in-sample LD reference? Otherwise consider using DENTIST!"))
} else {
  run_dentist(
    D=dataset_aligned,
    locus_chr=opt$chr,
    locus_start=opt$start,
    locus_end=opt$end,
    bfile=opt$bfile,
    maf.thresh=opt$maf,
    random.number=random.number,
    dentist.bin="/ssu/bsssu/software_bsssu/DENTIST"
  )
}


# Compute LD matrix
susie_ld <- prep_susie_ld(
  D=dataset_aligned,
  locus_chr=opt$chr,
  locus_start=opt$start,
  locus_end=opt$end,
  bfile=opt$bfile,
  maf.thresh=opt$maf,
  random.number=random.number,
  skip_dentist=opt$skip_dentist
)

# Compute trait variance
#dataset_aligned$MAF <- ifelse(dataset_aligned$freq < 0.5, dataset_aligned$freq, (1-dataset_aligned$freq))
D_var_y <- median(dataset_aligned$se^2*dataset_aligned$N*2*dataset_aligned$freq*(1-dataset_aligned$freq), na.rm = T)

# Filter full GWAS sum stat for locus region
D_sub <- dataset_aligned[match(rownames(susie_ld),dataset_aligned$SNP),]


### Run SUSIE

# 1) Check if susie output was produced
# 2) If susie output was produced, check that cs are not empty. If not, lower coverage until at least a cs is found or bottom threshold for coverage is reached
# 3) If cs are not empty, apply QC and then check that QCed cs object is not empty


coverage_value <- opt$cs_thresh
min_coverage <- 0.70
susie_error_message <- NULL  # Initialize an object to store the error message
fitted_rss <- NULL  # Initialize fitted_rss
max_iter=100*100

# While loop to adjust coverage until the condition is met
while (is.null(fitted_rss$sets$cs) && coverage_value >= min_coverage) {
  
  fitted_rss <- tryCatch( ### fitting susie with error catch
    {
      susie_rss(
        bhat = D_sub$b, 
        shat = D_sub$se, 
        n = max(D_sub$N), 
        R = susie_ld, 
        var_y = D_var_y, 
        L = 10,
        estimate_residual_variance = FALSE,
        coverage = coverage_value,
        max_iter = max_iter
      )
    },
    error = function(e) {
      susie_error_message <<- e$message  # Save the error message
      message("An error occurred: ", e$message)
      return(NULL)
    }
  )
  
  # Check for specific error message to exit - switch to cojo
  if (is.null(fitted_rss) && grepl("The estimated prior variance is unreasonably large", susie_error_message)) {
    write.table(susie_error_message, "failed_susie.txt", row.names = FALSE, col.names = FALSE)
    quit(save = "no", status = 0, runLast = FALSE)  # Exit the script gracefully
  }
  
  # Decrease coverage if fitted_rss$sets$cs is still null
  if (is.null(fitted_rss$sets$cs)) {
    coverage_value <- coverage_value - 0.01  # Decrease coverage by 0.01
    message("Trying again with coverage: ", coverage_value)
  }
}

# If still not cs after lowering coverage to 0.70, exit
if (is.null(fitted_rss$sets$cs)) {
  message(paste0("Final attempt failed, reached minimum coverage of ", min_coverage))
  quit(save = "no", status = 0, runLast = FALSE)  # Exit the script gracefully
} else {
  message("Success! Found sets with coverage: ", coverage_value)

  # Store final coverage and convergence status ???????????
  susie_final_coverage <- coverage_value 
  susie_convergence <- fitted_rss$converged ######################
  
#### Sodbo's function QCing Susie output --> check with him all the parameters required
  fitted_rss_cleaned <- susie.cs.ht(
    fitted_rss,
    D_sub$p,
    cs_lbf_thr = 2,
    signal_pval_threshold = 1e-6,
    purity_mean_r2_threshold = 0.5,
    purity_min_r2_threshold = 0.5,
    verbose = FALSE
  )
    
### Proceed only if fitted_rss_cleaned is not null    
  if(!is.null(fitted_rss_cleaned)){
    
    # vector of AF with SNP names
    freq <- D_sub$freq
    names(freq) <- D_sub$SNP
    
    # vector of sample size with SNP names
    N <- D_sub$N
    names(N) <- D_sub$SNP
    
    finemap.res <- lapply(fitted_rss_cleaned$sets$cs_index, function(x){
        
      # Extract lABF values  
      lABF_df <- data.frame(
        SNP = colnames(fitted_rss_cleaned$lbf_variable),
        lABF = fitted_rss_cleaned$lbf_variable[x,]
      )
      
      # Extract index of cs SNPs
      index <- paste0("L", x)
      cs_snps <- susie_get_cs(fitted_rss_cleaned, coverage=susie_final_coverage)$cs[[index]] ### INCLUDING COVERAGE IS CRUCIAL, OTHERWISE YOU GET DIFFERENT RESULTS!!!
      
      # Extract SNPs info, plus whether they're in the cs or not    
      susie_reformat <- D_sub
      susie_reformat$is_cs <- susie_reformat$SNP %in% colnames(fitted_rss_cleaned$lbf_variable)[cs_snps]
      
      susie_reformat <- susie_reformat %>%
        inner_join(lABF_df, by="SNP") %>%
        dplyr::select(SNP,BP,lABF,b,p,is_cs) %>%
        dplyr::rename(snp=SNP, position=BP) %>%
        dplyr::arrange(desc(lABF), p)
      
      effect = proccess_susie(fitted_rss_cleaned, x, freq, N)
      qc_metrics = fitted_rss_cleaned$sets$purity[paste0("L",x),]
        
      return(
        list(
          finemapping_lABFs = susie_reformat,
          effect = effect,
          qc_metrics = qc_metrics
          )
        )
    })
      
# Name set with the highest lABF SNP      
    names(finemap.res) <- sapply(finemap.res, function(x) x$finemapping_lABF$snp[1])
    
    
    #########################################
    # Organise list of what needs to be saved
    #########################################

    core_file_name <- paste0(opt$study_id, "_", ifelse(opt$phenotype_id == "full", gsub("_full", "", opt$phenotype_id), opt$phenotype_id))
    

    ## Create and save ind.snps-like table
    ind.snps <- D_sub %>%
      dplyr::filter(SNP %in% names(finemap.res)) %>%
      dplyr::mutate(freq_geno=NA, bJ=b, bJ_se=se, pJ=p, LD_r=NA, start=opt$start, end=opt$end, study_id=opt$study_id) %>%
      dplyr::select(CHR,SNP,BP,A1,freq,b,se,p,N,freq_geno,bJ,bJ_se,pJ,LD_r,snp_original,A2,type,any_of(c("sdY","s")),start,end,study_id,phenotype_id) %>%
      dplyr::rename("Chr"="CHR","bp"="BP","refA"="A1","n"="N","othA"="A2")
    
    fwrite(ind.snps, paste0(core_file_name, "_locus_chr", locus_name,"_ind_snps.tsv"), sep="\t", quote=F, na=NA)

    
    ## Save lABF of each conditional dataset
    lapply(names(finemap.res), function(x){
        
      sp_file_name <- paste0(core_file_name, "_", x, "_locus_chr", locus_name)
      
      # Create list object for lABFs to save in .rds file
      finemap_list <- finemap.res[[x]]
        
        # .rds object collecting 1) lABF, 2) pos for all SNPs, 3) list of SNPs in the credible set
      saveRDS(finemap_list, file = paste0(sp_file_name, "_susie_finemap.rds"))
        
        # .tsv with 1) study id and trait (if molQTL) locus info, 2) list of SNPs in the 99% credible set, 3) path and name of correspondent .rds file and 4) path and name of correspondent ind_snps.tsv table
        #  --> append each row to a master table collecting all info from processed sum stats
        ### Idea: create guidelines for generating study ids
        
      tmp <- data.frame(
        study_id = opt$study_id,
        phenotype_id = ifelse(opt$phenotype_id=="full", NA, opt$phenotype_id),
        credible_set = paste0(finemap_list$finemapping_lABFs %>% filter(is_cs==TRUE) %>% pull(snp), collapse=","),
        top_pvalue = min(finemap_list$finemapping_lABFs$p, na.rm=T),
          #### Nextflow working directory "work" hard coded - KEEP in mind!! #### 
        path_rds = paste0(opt$results_path, "/results/finemap/", sp_file_name, "_susie_finemap.rds"),
        path_ind_snps = paste0(opt$results_path, "/results/gwas_and_loci_tables/", opt$study_id, "_final_ind_snps_table.tsv"),
        chr=opt$chr
      )
      
      fwrite(tmp, paste0(sp_file_name, "_susie_coloc_info_table.tsv"), sep="\t", quote=F, col.names = F, na=NA)
    })
  }
}
