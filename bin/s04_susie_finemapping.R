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
  make_option("--susie_qc_cs_bf_thr", default=3, help="Credible set BF threshold for credible sets QC"),
  make_option("--susie_qc_pval_thr", default=1, help="Top SNP p-value threshold for credible sets QC"),
  make_option("--susie_qc_mean_r2_thr", default=0, help="Credible set purity mean r2 threshold for credible sets QC"),
  make_option("--susie_qc_min_r2_thr", default=0, help="Credible set for purity minimum r2 threshold for credible sets QC"),
  make_option("--publish_susie", default=FALSE, help=" Whether to publish the susie finemap .rds intermediate files"),
  make_option("--results_path", default=NULL, help="Path to \"/results\" folder"),
  make_option("--study_id", default=NULL, help="Id of the study"),
  make_option("--batch", default=NULL, help="File with multiple loci (chr, start, end, phenotype_id)")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

# Initialize a list to store all cleaned fitted_rss objects
all_fitted_rss_cleaned <- list()

# GWAS input
if (is.null(opt$batch)) { 
  opt$chr <- as.numeric(opt$chr)
  opt$start <- as.numeric(opt$start)
  opt$end <- as.numeric(opt$end)
  dataset_aligned <- fread(cmd=paste0("tabix ", opt$dataset_aligned, " ", opt$phenotype_id))
  colnames(dataset_aligned) <- c("phenotype_id", "snp_original","SNP","CHR","BP","A1","A2","freq","b","se","p","N", "type","temp")

  # Set up phenotypic variance correctly
  if(unique(dataset_aligned$type=="quant")){
    D_var_y = unique(dataset_aligned$temp)^2
  } else if(unique(dataset_aligned$type=="cc")){
    D_var_y = unique(dataset_aligned$temp) * (1-unique(dataset_aligned$temp))
  }
  dataset_aligned <- dataset_aligned |> dplyr::select(-temp)
  dataset_aligned[,"study_id"] = opt$study_id
  loci_df <- data.frame(
    chr = opt$chr,
    start = opt$start,
    end = opt$end,
    phenotype_id = opt$phenotype_id,
    study_id =  opt$study_id
    )
  } else {

  loci_df <- fread(opt$batch, header = TRUE)
  dataset_aligned <- fread(opt$dataset_aligned, header = TRUE)
  dataset_aligned <- dataset_aligned %>%
    rename(
      "p" = "SNP_PVAL",
      "freq" = "EAF",
      "b" = "BETA",
      "se" = "SE",
      "SNP"="SNPID",
      "BP" = "SNP_POS"
    )

  loci_df <- loci_df %>%
      rename(
        "chr" = "CHR",
        "start" = "START",
        "end" = "END"
      )
  if(loci_df[1,'TYPE'] == "gwas"){
      dataset_aligned <- dataset_aligned %>%
      rename(
        "study_id" = "TRAIT"
      )
      loci_df <- loci_df %>%
      rename(
        "study_id" = "TRAIT"
        )
      dataset_aligned = dataset_aligned %>% mutate("phenotype_id" = "full")
      loci_df = loci_df %>% mutate("phenotype_id" = "full")
      #loci_df["phenotype_id"]= "full"
      
    } else{
      loci_df <- loci_df %>%
        rename(
          "study_id" = "CELL",
          "phenotype_id" = "GENE"
        )
      dataset_aligned <- dataset_aligned %>%
      rename(
        "study_id" = "CELL",
        "phenotype_id" = "GENE"
      )
  }
  D_var_y = loci_df$PHENO_VAR
}

for (i in 1:nrow(loci_df)) {
  chr <- loci_df$chr[i]
  start <- loci_df$start[i]
  end <- loci_df$end[i]
  phenotype_name <- loci_df$phenotype_id[i]
  study_name <- loci_df$study_id[i] 
  
  locus_name <- paste0(chr, ":", start, ":", end)
  chr <- as.numeric(chr)
  start <- as.numeric(start)
  end <- as.numeric(end)
  dataset_aligned_sub <- dataset_aligned %>%
                filter((phenotype_id == phenotype_name) & (study_id == study_name))

  random.number <- stri_rand_strings(n=1, length=20, pattern="[A-Za-z0-9]")

### If required, run DENTIST to identify mismatches between GWAS sum stats and LD panel
  if (opt$skip_dentist){
    cat(paste0("I assume you provided in-sample LD reference? Otherwise consider using DENTIST!"))
  } else {
    run_dentist(
      D=dataset_aligned_sub,
      locus_chr=chr,
      locus_start=start,
      locus_end=end,
      bfile=opt$bfile,
      maf.thresh=opt$maf,
      random.number=random.number,
      dentist.bin="DENTIST"
    )
  }

  # Compute LD matrix
  susie_ld <- prep_susie_ld(
  D=dataset_aligned_sub,
  locus_chr=chr,
  locus_start=start,
  locus_end=end,
  bfile=opt$bfile,
  maf.thresh=opt$maf,
  random.number=random.number,
  skip_dentist=opt$skip_dentist
  )

# Check if susie_ld is NULL
  if (is.null(susie_ld)) {
    susie_error_message <- "prep_susie_ld returned NULL. No SNPs found in the locus or LD matrix could not be computed."
    
    no_variants_in_ld_ref <- data.frame(
      study_id = study_name,
      phenotype_id = phenotype_name,
      chr = chr,
      start = start,
      end = end,
      not_finemapped_reason = susie_error_message
    )

    fwrite(no_variants_in_ld_ref, paste0(random.number, "_NOT_FINEMAPPED_no_variants_from_locus_in_LD_ref.tsv"), sep="\t", na=NA, quote=F)
    next # skip this iteration and continue with the next one
  }

  # Filter full GWAS sum stat for locus region
  D_sub <- dataset_aligned_sub[match(rownames(susie_ld),dataset_aligned_sub$SNP),]


  ### Run SUSIE

  # 1) Check if susie output was produced
  # 2) If susie output was produced, check that cs are not empty. If not, lower coverage until at least a cs is found or bottom threshold for coverage is reached
  # 3) If cs are not empty, apply QC and then check that QCed cs object is not empty

  min_coverage <- 0.7
  L <- 10

  fitted_rss <- run_susie_w_retries(
    D_sub,
    D_var_y[i],
    susie_ld,
    L = L,
    coverage = opt$cs_thresh,
    min_coverage = min_coverage,
    max_iter = opt$susie_max_iter,
    min_abs_corr = NULL
  )

# If successfull, pass to QC
if (!is.null(fitted_rss) && !is.null(fitted_rss$sets$cs)) {
  
  # Perform QC only for loci fine-mapped with L=10 (or whatever number was set) or, if N SNPs < L, with L=N SNPs in the locus
  if (length(fitted_rss$KL) == L | length(fitted_rss$KL) == nrow(D_sub)){

  ### Post susie QC od credible sets
  #  fitted_rss_cleaned <- flanders::susie.cs.ht( ### THIS IS TEMPORARY UNTIL UPDATE OF THE GITHUB FLANDERS R REPO
    fitted_rss_cleaned <- susie.cs.ht(
      fitted_rss,
      D_sub$p,
      cs_bf_thr = opt$susie_qc_cs_bf_thr,
      signal_pval_threshold = opt$susie_qc_pval_thr,
      purity_mean_r2_threshold = opt$susie_qc_mean_r2_thr,
      purity_min_r2_threshold = opt$susie_qc_min_r2_thr,
      verbose = TRUE
    )
    fitted_rss_cleaned$metadata <- loci_df[i, ] |> select(study_id, phenotype_id, chr, start, end, TYPE, N)

    ### Re-finemap with L=1 if all cs gets removed by QC
      if(is.null(fitted_rss_cleaned)){
        fitted_rss_cleaned <- run_susie_w_tryCatch( ### Make sure this works as intended
          D_sub,
          D_var_y[i],
          susie_ld,
          L = 1,
          coverage = opt$cs_thresh,
          max_iter = opt$susie_max_iter,
          min_abs_corr = 0
        )
        fitted_rss_cleaned$comment_section <- paste0("Locus re-finemapped at L=1 after none of the credible sets fine-mapped at L=", L, " passed post susie QC")
        fitted_rss_cleaned$metadata <- loci_df[i, ] |> select(study_id, phenotype_id, chr, start, end, TYPE, N)
      }

    # Skip QC for loci fine-mapped with L=1
    } else if (length(fitted_rss$KL)==1){
      fitted_rss_cleaned <- fitted_rss
      fitted_rss_cleaned$metadata <- loci_df[i, ] |> select(study_id, phenotype_id, chr, start, end, TYPE, N)
    }

    ## List of loci which were still fine-mapped but with L=1 (and why)
    if(!is.na(fitted_rss_cleaned$comment_section)){
      L1_finemap <- data.frame(
        study_id = study_name,
        phenotype_id = phenotype_name,
        chr = chr,
        start = start,
        end = end,
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

      L1_finemap_missing_loci_post_susie_QC <- L1_finemap |> dplyr::filter(grepl(paste0("Locus re-finemapped at L=1 after none of the credible sets fine-mapped at L=", L," passed post susie QC"), finemapped_L1_reason))
      if(nrow(L1_finemap_missing_loci_post_susie_QC) > 0){
        fwrite(L1_finemap_missing_loci_post_susie_QC, paste0(random.number, "_FINEMAPPED_L1_recover_after_susie_QC.tsv"), sep="\t", na=NA, quote=F)
      }
    
    }

  all_fitted_rss_cleaned[[locus_name]] <- fitted_rss_cleaned ### store all in a list

  } else { ### if region was not fine-mapped at all!
    
    failed_finemap <- data.frame(
      study_id = study_name,
      phenotype_id = phenotype_name,
      chr = chr,
      start = start,
      end = end
    )
    fwrite(failed_finemap, paste0(random.number, "_NOT_FINEMAPPED_no_credible_sets_found.tsv"), sep="\t", na=NA, quote=F)
  }
}

# Expand credible sets
expanded_cs <- lapply(all_fitted_rss_cleaned, expand_cs)
names(expanded_cs) <- names(all_fitted_rss_cleaned)

# Add conditional beta and se 
beta_se_list <- lapply(all_fitted_rss_cleaned, function(locus){
  lapply(locus$sets$cs_index, function(x){
    get_beta_se_susie(locus, x)
  })
})

### From Susie object to annData
ad <- from_susie_to_anndata(
  finemap_list = all_fitted_rss_cleaned, 
  cs_indices = expanded_cs,
  beta_se_cond = beta_se_list
)

message("AnnData created!")

## Save anndata
anndata::write_h5ad(ad, paste0("batch_", opt$batch, "_anndata.h5ad"))

