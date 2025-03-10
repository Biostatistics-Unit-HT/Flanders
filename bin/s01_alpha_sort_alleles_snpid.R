suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
    make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
    make_option("--bfile", default=NULL, help="Path and prefix name of custom LD bfiles (PLINK format .bed .bim .fam)"),
    make_option("--grch", default=NULL, help="Genome reference build of GWAS sum stats")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

## Change snpids
bim <- fread(paste0(opt$bfile, ".bim"))
names(bim) <- c("CHR","snp_original","V3", "BP","V5","V6")

# If necessary, lift to build 38
if(as.numeric(opt$grch)==37){
  
  bim_lifted <- hg19ToHg38_liftover(bim)
  
# Remove rows with duplicated SNP (all occurrences!)
  bim_lifted_no_dups <- bim_lifted[, if (.N == 1) .SD, by = snp_original]

# Save list of SNP ids to extract from .bed
  fwrite(bim_lifted_no_dups %>% dplyr::select(snp_original), paste0(opt$bfile, "_snps_to_extract.txt"), col.names=F, quote=F)
  
# Extract list of SNPs from .bim (to match it with .bed!)
  system(paste0("plink2 --bfile ", opt$bfile, " --extract ", opt$bfile, "_snps_to_extract.txt --make-bed --out ", opt$bfile, "_alpha_sorted_alleles"))
  system(paste0("rm ", opt$bfile, "_snps_to_extract.txt"))
  
  bim_cleaned <- bim_lifted_no_dups
  
} else if(as.numeric(opt$grch)==38){

  # Since you're dong it for grch 37 bfiles, do it also for grch 38 - create sym links of .bed and .fam
  bim_cleaned <- bim
  system(paste0("ln -s ", opt$bfile, ".fam ", opt$bfile, "_alpha_sorted_alleles.fam"))
  system(paste0("ln -s ", opt$bfile, ".bed ", opt$bfile, "_alpha_sorted_alleles.bed"))
}

# Alpha sort alleles
bim_alpha_sorted <- bim_cleaned %>%
  dplyr::mutate(
    A1 = pmin(V5, V6), # Sort A1 and A2 alphabetically
    A2 = pmax(V5, V6),
    snp_original = paste0("chr", CHR, ":", BP, ":", A1, ":", A2)
  ) %>%
  dplyr::select(CHR, snp_original, V3, BP, V5, V6)

fwrite(
  bim_alpha_sorted,
  paste0(gsub(".*/(.*).bim", "\\1", opt$bfile), "_alpha_sorted_alleles.bim"),
  quote=F, na=NA, sep="\t", col.names = F)

