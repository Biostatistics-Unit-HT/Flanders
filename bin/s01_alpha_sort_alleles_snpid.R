#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
liftOver <- rtracklayer::liftOver
import.chain <- rtracklayer::import.chain
GRanges <- GenomicRanges::GRanges
IRanges <- IRanges::IRanges

hg19ToHg38_liftover <- function(
    dataset_munged,
    default_chain_file = "hg19ToHg38.over.chain"
){
  
  if(!(file.exists(default_chain_file))){
    ### Download chain file and unzip it
    exit_status = system("wget -O - http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz | gunzip -c > hg19ToHg38.over.chain")

    # Raise an error if the external command fails
    if (exit_status != 0) {
      cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
      quit(status = 1, save = "no")
    }
  }

  ch <- import.chain(default_chain_file)

  dt_for_ranges <- copy(dataset_munged)
  dt_for_ranges[, start := BP]
  dt_for_ranges <- unique(dt_for_ranges, by = c("snp_original", "CHR", "start", "BP"))
  
  dataset_ranges <- GRanges(
    seqnames = paste0("chr", dt_for_ranges$CHR),
    ranges = IRanges(start = dt_for_ranges$start, end = dt_for_ranges$BP),
    snp_original = dt_for_ranges$snp_original
  )
  
  rm(dt_for_ranges)
  gc()
  
  dataset_ranges38 <- liftOver(dataset_ranges, ch)
  dataset_ranges38_df <- as.data.table(unlist(dataset_ranges38))
  setnames(dataset_ranges38_df, "end", "BP")
  dataset_ranges38_df <- dataset_ranges38_df[, .(BP, snp_original)]
  
  dataset_lifted <- merge(dataset_munged[, !"BP"], dataset_ranges38_df, by = "snp_original", all = FALSE)
  
  return(dataset_lifted)
}

# Get arguments specified in the sbatch
option_list <- list(
    make_option("--bfile", default=NULL, help="Path and prefix name of custom LD bfiles (PLINK format .bed .bim .fam)"),
    make_option("--run_liftover", type = "logical", default=TRUE, help="Perform liftover to GRCh38?"),
    make_option("--grch", default=NULL, help="Genome reference build of GWAS sum stats")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Change snpids
bim <- fread(paste0(opt$bfile, ".bim"))
names(bim) <- c("CHR","snp_original","V3", "BP","V5","V6")

# If necessary, lift to build 38
if(as.numeric(opt$grch)==37 && as.logical(opt$run_liftover)){
  
  bim_lifted <- hg19ToHg38_liftover(bim)
  
# Remove rows with duplicated SNP (all occurrences!)
  bim_lifted_no_dups <- bim_lifted[, if (.N == 1) .SD, by = snp_original]

# Save list of SNP ids to extract from .bed
  fwrite(bim_lifted_no_dups %>% dplyr::select(snp_original), paste0(opt$bfile, "_snps_to_extract.txt"), col.names=F, quote=F)
  
# Extract list of SNPs from .bim (to match it with .bed!)
  exit_status = system(paste0("plink2 --bfile ", opt$bfile, " --extract ", opt$bfile, "_snps_to_extract.txt --make-bed --out ", opt$bfile, ".GRCh38"))
  
  # Raise an error if the external command fails
  if (exit_status != 0) {
    cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
    quit(status = 1, save = "no")
  }
  
  system(paste0("rm ", opt$bfile, "_snps_to_extract.txt"))
  
  bim_cleaned <- bim_lifted_no_dups
  
} else if(as.numeric(opt$grch)==38){
  
  bim_cleaned <- bim
  system(paste0("cp ", opt$bfile, ".fam ", opt$bfile, ".GRCh38.alpha_sorted_alleles.fam"))
  system(paste0("cp ", opt$bfile, ".bed ", opt$bfile, ".GRCh38.alpha_sorted_alleles.bed"))
  #system(paste0("ln -s ", opt$bfile, ".fam ", opt$bfile, "_alpha_sorted_alleles.fam"))
  #system(paste0("ln -s ", opt$bfile, ".bed ", opt$bfile, "_alpha_sorted_alleles.bed"))
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
  paste0(basename(opt$bfile), ".GRCh38.alpha_sorted_alleles.bim"),
  quote=F, na=NA, sep="\t", col.names = F)

