#!/usr/bin/env -S Rscript --vanilla

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
  make_option("--susie_max_iter", default=400, help="Maximum number of susie iterations"),
  make_option("--publish_susie", default=FALSE, help=" Whether to publish the susie finemap .rds intermediate files"),
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

# Set up phenotypic variance correctly
if(unique(dataset_aligned$type=="quant")){
  D_var_y = unique(dataset_aligned$temp)^2
} else if(unique(dataset_aligned$type=="cc")){
  D_var_y = unique(dataset_aligned$temp) * (1-unique(dataset_aligned$temp))
}
dataset_aligned <- dataset_aligned |> dplyr::select(-temp)

################################
# Conditional analysis with SUSIE
################################

random.number <- stri_rand_strings(n=1, length=20, pattern="[A-Za-z0-9]")

### If required, run DENTIST to identify mismatches between GWAS sum stats and LD panel
if (opt$skip_dentist){
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
    dentist.bin="DENTIST"
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

# Filter full GWAS sum stat for locus region
D_sub <- dataset_aligned[match(rownames(susie_ld),dataset_aligned$SNP),]


### Run SUSIE

# 1) Check if susie output was produced
# 2) If susie output was produced, check that cs are not empty. If not, lower coverage until at least a cs is found or bottom threshold for coverage is reached
# 3) If cs are not empty, apply QC and then check that QCed cs object is not empty

min_coverage <- 0.7
L <- 10

fitted_rss <- run_susie_w_retries(
  D_sub,
  D_var_y,
  susie_ld,
  L = L,
  coverage = opt$cs_thresh,
  min_coverage = min_coverage,
  max_iter = opt$susie_max_iter,
  min_abs_corr = NULL
)

# If successful, perform QC
if (!is.null(fitted_rss) && !is.null(fitted_rss$sets$cs)) {

#### Sodbo's function QCing Susie output --> check with him all the parameters required
  fitted_rss_cleaned <- susie.cs.ht(
    fitted_rss,
    D_sub$p,
    cs_lbf_thr = 2, # TODO: this should be a pipeline argument
    signal_pval_threshold = 1, # TODO: this should be a pipeline argument
    purity_mean_r2_threshold = 0.5, # TODO: this should be a pipeline argument
    purity_min_r2_threshold = 0.5, # TODO: this should be a pipeline argument
    verbose = TRUE
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
      
      beta_se_list <- get_beta_se_susie(fitted_rss_cleaned, x)
        
      # Extract lABF values  
      lABF_df <- data.frame(
        SNP = colnames(fitted_rss_cleaned$lbf_variable),
        lABF = fitted_rss_cleaned$lbf_variable[x,],
        bC = beta_se_list$beta,
        bC_se = beta_se_list$se
      )
      
      # Extract index of cs SNPs
      index <- paste0("L", x)
      cs_snps <- susie_get_cs(fitted_rss_cleaned, coverage=fitted_rss_cleaned$sets$requested_coverage)$cs[[index]] ### INCLUDING COVERAGE IS CRUCIAL, OTHERWISE YOU GET DIFFERENT RESULTS!!!
      
      # Extract SNPs info, plus whether they're in the cs or not    
      susie_reformat <- D_sub
      susie_reformat$is_cs <- susie_reformat$SNP %in% colnames(fitted_rss_cleaned$lbf_variable)[cs_snps]
      
      susie_reformat <- susie_reformat |>
        dplyr::inner_join(lABF_df, by="SNP") |>
        dplyr::select(SNP,BP,lABF,bC,bC_se,is_cs) |>
        dplyr::rename(snp=SNP, position=BP) |>
        dplyr::arrange(desc(lABF))
      
      effect <- susie_reformat[1,] # select the first row which represents the
      # top SNP by lABF. 
      snp_top <- effect$snp
      chr_pos_a1_a0_top <- strsplit(effect$snp, ":")[[1]]
      a1_top <- chr_pos_a1_a0_top[3]
      a0_top <- chr_pos_a1_a0_top[4]
      freq_top <- freq[snp_top]
      N_top <- N[snp_top]
      beta_top <- effect$bC
      se_top <- effect$bC_se
      
      effect <- data.frame(
        snp = snp_top,
        a1 = a1_top,
        a0 = a0_top,
        freq = freq_top,
        N = N_top,
        beta = beta_top,
        se = se_top
        )
      
      qc_metrics <- fitted_rss_cleaned$sets$purity[paste0("L",x),] |>
        dplyr::mutate(coverage = fitted_rss_cleaned$sets$coverage[x], L = length(fitted_rss_cleaned$KL)) # Add also requested coverage and L
      
      metadata_df <- data.frame(
        study_id=opt$study_id,
        phenotype_id=opt$phenotype_id,
        chr=opt$chr,
        start=opt$start,
        end=opt$end
      )
      
      return(
        list(
          finemapping_lABFs = susie_reformat,
          effect = effect,
          qc_metrics = qc_metrics,
          metadata = metadata_df
          )
        )
    })

    # Name each credible set
    names(finemap.res) <- paste(
      paste0("chr", opt$chr),
      opt$study_id,
      opt$phenotype_id,
      sapply(finemap.res, function(x) x$finemapping_lABF$snp[1]),
      sep="::"
    )
    
    #########################################
    # Organise list of what needs to be saved
    #########################################

    core_file_name <- paste0(opt$study_id, "_", opt$phenotype_id)
#    if(opt$phenotype_id=="full") { core_file_name <- gsub("_full", "", core_file_name)}

    ## Save .rds object
    saveRDS(finemap.res, file = paste0(core_file_name, "_locus_chr", locus_name, "_susie_finemap.rds"))

    ## Save info about each cs
    tmp <- rbindlist(lapply(finemap.res, function(x){              
      data.frame(
        credible_set_name = paste(
          paste0("chr", opt$chr),
          opt$study_id,
          opt$phenotype_id,
          x$finemapping_lABF$snp[1],
          sep="::"
        ),
        credible_set_snps = paste0(x$finemapping_lABFs |> dplyr::filter(is_cs==TRUE) |> dplyr::pull(snp), collapse=","),
        study_id = opt$study_id,
        phenotype_id = opt$phenotype_id,
        chr = opt$chr,
        start = opt$start,
        end = opt$end,
        top_pvalue = min(pchisq((x$finemapping_lABFs$bC/x$finemapping_lABFs$bC_se)**2, 1, lower.tail=FALSE), na.rm=T),
          #### Nextflow working directory "work" hard coded - KEEP in mind!! #### 
        path_rds = ifelse(
          opt$publish_susie,
          paste0(opt$results_path, "/results/finemap/", core_file_name, "_locus_chr", locus_name, "_susie_finemap.rds"),
          NA),
        x$effect
      ) |>
        dplyr::rename(bC=beta, bC_se=se)
    }))
    fwrite(tmp, paste0(core_file_name, "_locus_chr", locus_name, "_cs_info_table.tsv"), sep="\t", quote=F, col.names = F, na=NA)
    
    
    ## List of loci which were still fine-mapped but with L=1 (and why)
    if(!is.na(fitted_rss_cleaned$comment_section)){
      L1_finemap <- data.frame(
        study_id = opt$study_id,
        phenotype_id = opt$phenotype_id,
        chr = opt$chr,
        start = opt$start,
        end = opt$end,
        finemapped_L1_reason = fitted_rss_cleaned$comment_section
      )
      
      L1_finemap_variance_too_large <- L1_finemap |> dplyr::filter(grepl("The estimated prior variance is unreasonably large", finemapped_L1_reason))
      if(nrow(L1_finemap_variance_too_large) > 0){
        fwrite(L1_finemap_variance_too_large, paste0(random.number, "_FINEMAPPED_L1_prior_variance_too_large.tsv"), sep="\t", na=NA, quote=F)
      }
      
      L1_finemap_did_not_converge <- L1_finemap |> dplyr::filter(grepl("IBSS algorithm did not converge", finemapped_L1_reason))
      if(nrow(L1_finemap_did_not_converge) > 0){
        fwrite(L1_finemap_did_not_converge, paste0(random.number, "_FINEMAPPED_L1_IBSS_algorithm_did_not_converge.tsv"), sep="\t", na=NA, quote=F)
      }      
  
    }
    
  }
} else { ### if region was not fine-mapped at all!
  
  failed_finemap <- data.frame(
    study_id = opt$study_id,
    phenotype_id = opt$phenotype_id,
    chr = opt$chr,
    start = opt$start,
    end = opt$end
  )
  fwrite(failed_finemap, paste0(random.number, "_NOT_FINEMAPPED_no_credible_sets_found.tsv"), sep="\t", na=NA, quote=F)
  
}
