#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))

option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--dataset_aligned", default=NULL, help="GWAS sum stats munged and aligned"),
  make_option("--maf", default=1e-04, help="MAF filter", metavar="character"),
  make_option("--p_thresh1", default=5e-08, help="Significant p-value threshold for top hits"),
  make_option("--p_thresh2", default=1e-05, help="P-value threshold for loci borders"),  
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--study_id", default=NULL, help="Id of the study")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

# Load-in summary statistics munged and aligned and filter by MAF (doesn't make sense to do this later!)
dataset_aligned <- fread(opt$dataset_aligned, data.table=F) %>%
  dplyr::filter(MAF > opt$maf) %>%
  dplyr::arrange(CHR)


################
# Locus breaker
################

loci_list <- as.data.frame(rbindlist(
  lapply(dataset_aligned %>% group_split(phenotype_id), function(x){
    ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
    if(any(x %>% pull(p) < opt$p_thresh1)){
      ### Loci identification
      locus.breaker(
        x,
        p.sig=opt$p_thresh1,
        p.limit=opt$p_thresh2,
        hole.size=opt$hole,
        p.label="p",
        chr.label="CHR",
        pos.label="BP")
    }
  }) %>% discard(is.null)
))

# Slightly enlarge locus by 200kb!
loci_list$start <- as.numeric(loci_list$start) -100000
loci_list$start[loci_list$start < 0] <- 0 # Replace negative numbers with 0
loci_list$end <- as.numeric(loci_list$end) + 100000


if(nrow(loci_list) > 0){
### Add study ID to the loci table. Save
  loci_list <- loci_list %>% mutate(study_id=opt$study_id)
  
  ### Check if locus spans the HLA locus chr6:28,510,120-33,480,577
  ### https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
  hla_start=28510120
  hla_end=33480577
  hla_coord <- seq(hla_start,hla_end)
  
  loci_list <- loci_list %>%
    mutate(is_in_hla = chr == 6 & ((start <= hla_end & end >= hla_start) | (start >= hla_start & end <= hla_end) | (start <= hla_end & end >= hla_end) | (start <= hla_start & end >= hla_start)))
  
  # This code checks for four conditions to determine if a locus partially or completely spans the HLA region:
  # 1) Locus starts before (or at the exact beginning of) HLA and ends after (or at the exact end of) HLA.
  # 2) Locus starts and ends within the HLA region.
  # 3) Locus starts before the end of HLA and ends after the end of HLA.
  # 4) Locus starts before the beginning of HLA and ends after the beginning of HLA.
  
  cat(paste0("\n", nrow(loci_list), " significant loci identified for ", opt$study_id, "\n"))
  fwrite(loci_list, paste0(opt$study_id, "_loci.tsv"), sep="\t", quote=F, na=NA)
}



