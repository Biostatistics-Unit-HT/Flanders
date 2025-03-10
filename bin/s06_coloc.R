suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--coloc_guide_table", default=NULL, help="Path and filename of table listing all coloc pairwise tests to perform"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis"),
  make_option("--chr_cs", default=NULL, help="Chromosomes to split the analysis in")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))



###### COLOCALISTION ######

# Import table listing all coloc pairwise tests to perform
coloc_combo <- fread(opt$coloc_guide_table, data.table = F)

# Split the dataframe into a list of rows
coloc_combo_ls <- split(coloc_combo, seq(nrow(coloc_combo)))

# Perform coloc! 
coloc.full <- lapply(coloc_combo_ls, function(x){
  
# Load-in precomputed lABF  
  conditional.dataset1 <- readRDS(x$t1_path_rds)$finemapping_lABFs
  conditional.dataset2 <- readRDS(x$t2_path_rds)$finemapping_lABFs
  
# Retrieve important info from file name
  t1 <- ifelse(is.na(x$t1_phenotype_id), x$t1_study_id, paste0(x$t1_study_id, "_", x$t1_phenotype_id))
  top_snp1 <- gsub(paste0(".*/", t1, "_(.*)_locus_.*_finemap.rds"), "\\1", x$t1_path_rds)
  
  t2 <- ifelse(is.na(x$t2_phenotype_id), x$t2_study_id, paste0(x$t2_study_id, "_", x$t2_phenotype_id))
  top_snp2 <- gsub(paste0(".*/", t2, "_(.*)_locus_.*_finemap.rds"), "\\1", x$t2_path_rds)
  
# Perform colocalisation for each combination of independent SNPs
  coloc.res <- hcolo.cojo.ht(
    df1 = conditional.dataset1 %>% dplyr::select(snp, lABF),
    df2 = conditional.dataset2 %>% dplyr::select(snp, lABF)
  )
  
  # Add top SNPs and traits
  coloc.res$summary <- coloc.res$summary %>%
    mutate(
      t1_study_id=x$t1_study_id,
      t1=t1,
      hit1=top_snp1,
      t2_study_id=x$t2_study_id,
      t2=t2,
      hit2=top_snp2
    )
  
  coloc.res
})

# Store ALL the summary output in a data frame, adding tested traits column and SAVE 
only_summary_df <- as.data.frame(rbindlist(lapply(coloc.full, function(x) { x$summary })))
    
# Give random index for avoiding overwrite (?) --> SHOULD NOT BE NEEDED, TAKEN CARE BY NEXTFLOW
#random.number=stri_rand_strings(n=1, length=20, pattern = "[A-Za-z0-9]")
fwrite(only_summary_df, paste0(#random.number, "_", 
  opt$coloc_id, "_chr", opt$chr_cs,"_colocalization.table.all.tsv"), quote=F, sep="\t", na=NA)
