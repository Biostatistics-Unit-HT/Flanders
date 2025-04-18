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
  make_option("--p_thresh3", default=1e-04, help="Noise p-values threshold for COJO"),
  make_option("--maf", default=1e-04, help="MAF filter", metavar="character"),
  make_option("--bfile", default=NULL, help="Path and prefix name of custom LD bfiles (PLINK format .bed .bim .fam)"),
  make_option("--skip_dentist", default=TRUE, help="Whether to skip the match of SNPs LD between GWAS sum stat and LD reference (performed by DENTIST), and consequent removal of mismatched SNPs"),
  make_option("--p_thresh4", default=1e-06, help="P-value significant threshold for redefining loci boundaries post-COJO"),  
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
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
  dataset_aligned <- dataset_aligned %>% rename(sdY=temp)
} else {
  dataset_aligned <- dataset_aligned %>% rename(s=temp)
}


################################
# Conditional analysis with COJO
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
    dentist.bin="DENTIST"
  )
}


### Conditional analysis by 2-steps COJO - If DENTIST was performed, input file are already ready. If not, to be created

conditional.dataset <- cojo.ht(
  D=dataset_aligned,
  locus_chr=opt$chr,
  locus_start=opt$start,
  locus_end=opt$end,
  p.thresh=as.numeric(opt$p_thresh3),
  maf.thresh=as.numeric(opt$maf),
  bfile=opt$bfile,
  random.number=random.number,
  skip_dentist=opt$skip_dentist
)


##### if conditional.dataset is null (no independently associated SNP found), stop here

if(!is.null(conditional.dataset)){ ### is there smarter way to exit the script here, without throwing an error?

  # Plot conditioned GWAS sum stats
  pdf(paste0(opt$study_id, "_locus_chr", locus_name, "_conditioned_loci.pdf"), height=3.5*nrow(conditional.dataset$ind.snps), width=10) ### have the original loci boundaries in the name, or the slightly enlarged ones?
  plot.cojo.ht(conditional.dataset) + plot_annotation(paste("Locus chr", locus_name))
  dev.off()
  
  
  
  ####################
  # Locus breaker BIS
  ###################
  
  # Add (original) start and end columns to independent SNPs df
  conditional.dataset$ind.snps <- conditional.dataset$ind.snps %>%
    dplyr::mutate(
      start = opt$start,
      end = opt$end,
      study_id = opt$study_id,
      phenotype_id = opt$phenotype_id)
  
  # Recompute loci boundaries
  for (i in 1:length(conditional.dataset$results)) {
  
  ### Check if there's any SNP at p-value lower than the set threshold
    x <- conditional.dataset$results[[i]]
    
    if(isTRUE(any(x %>% pull(pC) < opt$p_thresh4))){
      new_bounds <- locus.breaker(
        x,
        p.sig = as.numeric(opt$p_thresh4),
        p.limit = as.numeric(opt$p_thresh3),
        hole.size = opt$hole,
        p.label="pC",
        chr.label="Chr",
        pos.label="bp")
      
      # Slightly enlarge locus by 200kb!
      new_bounds <- new_bounds %>%
        dplyr::mutate(start=as.numeric(start)-100000, end=as.numeric(end)+100000) %>%
        # It could happen (should be rare) that more than one locus if found - take only the one with lowest pC
        dplyr::arrange(pC) %>%
        slice(1)
      
      # Remove SNPs not included in loci boundaries
      conditional.dataset$results[[i]] <- x %>% 
        dplyr::filter(bp >= new_bounds$start & bp <= new_bounds$end)
      
      # Add locus start and end to ind snps df
      conditional.dataset$ind.snps[i,] <- conditional.dataset$ind.snps[i,] %>%
        dplyr::mutate(start=new_bounds$start, end=new_bounds$end)
    }
  }
  
  
  # Prepare table of independent SNPs - add study and pheno ids, whether a credible set will be calculated for each SNP
  #conditional.dataset$ind.snps <- conditional.dataset$ind.snps %>%
  #  dplyr::mutate(
      # Check whether the corresponding conditional dataset is empty (thus not passing the p-value thresholds)
  #    is_cs_avail = !(sapply(conditional.dataset$results, is.null))
  #  )
  
  # Remove eventually empty dataframes
  conditional.dataset$results <- conditional.dataset$results %>% discard(is.null)
  
  
  
  #############
  # Finemapping
  #############
  
  # Perform finemapping of each conditional dataset
  finemap.res <- lapply(conditional.dataset$results, function(x){
    finemap.cojo(x, cs_threshold=opt$cs_thresh)
  })
  
  
  
  #########################################
  # Organise list of what needs to be saved
  #########################################
  
  ## Save independent association signals
  core_file_name <- paste0(opt$study_id, "_", opt$phenotype_id)
  if(opt$phenotype_id=="full") { core_file_name <- gsub("_full", "", core_file_name)}
  
  fwrite(conditional.dataset$ind.snps, paste0(core_file_name, "_locus_chr", locus_name,"_ind_snps.tsv"), sep="\t", quote=F, na=NA)
  
  
  ## Save lABF of each conditional dataset
  lapply(finemap.res, function(x){
    
    sp_file_name <- paste0(core_file_name, "_", unique(x$cojo_snp), "_locus_chr", locus_name)
    
    # .rds object collecting 1) lABF, 2) beta, 3) pos for all SNPs, 3) list of SNPs in the credible set
    saveRDS(list(finemapping_lABFs=x) #%>% select(-cojo_snp)
            , file=paste0(sp_file_name, "_cojo_finemap.rds")) ### cojo_snp reported in the file name
    
    # .tsv with 1) study id and trait (if molQTL) locus info, 2) list of SNPs in the 99% credible set, 3) path and name of correspondent .rds file and 4) path and name of correspondent ind_snps.tsv table
    #  --> append each row to a master table collecting all info from processed sum stats
    ### Idea: create guidelines for generating study ids
    
    tmp <- data.frame(
      study_id = opt$study_id,
      phenotype_id = ifelse(opt$phenotype_id=="full", NA, opt$phenotype_id),
      credible_set = paste0(x %>% filter(is_cs==TRUE) %>% pull(snp), collapse=","),
      top_pvalue = min(x$pC, na.rm=T),
      path_rds = paste0(opt$results_path, "/results/finemap/", sp_file_name, "_cojo_finemap.rds"),
      path_ind_snps = paste0(opt$results_path, "/results/gwas_and_loci_tables/", opt$study_id, "_final_ind_snps_table.tsv"),
      chr=opt$chr
    )
    fwrite(tmp, paste0(sp_file_name, "_cojo_coloc_info_table.tsv"),
           sep="\t", quote=F, col.names = F, na=NA)
  })
}