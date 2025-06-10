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
bim <- readr::read_delim(paste0(opt$bfile, ".bim"))
bim <- as.data.table(bim)
names(bim) <- c("CHR","snp_original","V3", "BP","V5","V6")

# Make a standardized snp id as CHR:BP:V5:V6 and save this as reference for downstream operations
bim <- bim |>
  dplyr::mutate(
    snp_original = paste0(CHR, ":", BP, ":", V5, ":", V6)
  )
fwrite(
  bim,
  paste0(opt$bfile, ".standard_snpid.bim"),
  quote=F, na=NA, sep="\t", col.names = F
)

# Remove SNPs where SNP ID is duplicated since this can mess up liftOver and other operations downstream
if (sum(duplicated(bim$snp_original)) > 0 ) {
  message(sum(duplicated(bim$snp_original)), " duplicated snp ids in bim file. Removing them")
  bim <- bim[!(duplicated(bim[, .(snp_original)]) | duplicated(bim[, .(snp_original)], fromLast=T))]
}

# If necessary, lift to build 38
if(as.numeric(opt$grch)==37 && as.logical(opt$run_liftover)){
  
  bim_to_clean <- hg19ToHg38_liftover(bim)
  
} else if(as.numeric(opt$grch)==38){

  bim_to_clean <- bim

}
# Remove rows with duplicated SNP by CHR POS
# This get rid of multi-allelic variants and any other odd situations
bim_cleaned <- bim_to_clean[!(duplicated(bim_to_clean[, .(CHR, BP)]) | duplicated(bim_to_clean[, .(CHR, BP)], fromLast=T))]

# Save list of SNP ids to extract from .bed
extract_file <- paste0(opt$bfile, "_snps_to_extract.txt")
fwrite(list(bim_cleaned |> dplyr::pull(snp_original) |> unique()), extract_file, col.names=F, quote=F)

# Extract list of SNPs from .bim (to match it with .bed!)
exit_status = system(paste0("plink2 --bed ", opt$bfile, ".bed --fam ", opt$bfile, ".fam --bim ", opt$bfile, ".standard_snpid.bim --extract ", extract_file, " --make-bed --out ", opt$bfile, ".GRCh38.alpha_sorted_alleles"))
  
# Raise an error if the external command fails
if (exit_status != 0) {
    cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
    quit(status = 1, save = "no")
  }
  
# Alpha sort alleles
bim_alpha_sorted <- bim_cleaned |>
  dplyr::mutate(
    A1 = pmin(V5, V6), # Sort A1 and A2 alphabetically
    A2 = pmax(V5, V6),
    snp_original = paste0("chr", CHR, ":", BP, ":", A1, ":", A2)
  ) |>
  dplyr::select(CHR, snp_original, V3, BP, V5, V6)

fwrite(
  bim_alpha_sorted,
  paste0(opt$bfile, ".GRCh38.alpha_sorted_alleles.bim"),
  quote=F, na=NA, sep="\t", col.names = F)

