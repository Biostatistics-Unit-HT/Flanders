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
  
#  dataset_lifted <- merge(dataset_munged[, !"BP"], dataset_ranges38_df, by = "snp_original", all = FALSE)

# New merge which keeps original file order - CRUCIAL TO MANTAIN ALIGNMENT WITH .BED!!
  dataset_lifted <- dataset_ranges38_df[
    dataset_munged[, !"BP"],
    on = "snp_original",
    nomatch = 0
  ]

  return(dataset_lifted)
}

# Get arguments specified in the sbatch
option_list <- list(
    make_option("--bfile", default=NULL, help="Path and prefix name of custom LD genotype files (PLINK bed/bim/fam or pgen/psam/pvar)"),
    make_option("--file-type", default = NULL, help = "Input file type: plink1 or plink2"),
    make_option("--run_liftover", type = "logical", default=TRUE, help="Perform liftover to GRCh38?"),
    make_option("--grch", default=NULL, help="Genome reference build of GWAS sum stats")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Detect input format
is_pgen <- file.exists(paste0(opt$bfile, ".pgen"))

## Read variant info
if (is_pgen) {
  # .pvar has header lines starting with ## and a #CHROM header line
  pvar <- fread(paste0(opt$bfile, ".pvar"), skip = "#CHROM")
  names(pvar)[1] <- "CHR"
  bim <- data.table(
    CHR = pvar$CHR,
    snp_original = pvar$ID,
    V3 = 0,
    BP = pvar$POS,
    V5 = pvar$REF,
    V6 = pvar$ALT
  )
} else {
  bim <- fread(paste0(opt$bfile, ".bim"))
  names(bim) <- c("CHR","snp_original","V3", "BP","V5","V6")
}

# Make a standardized snp id as CHR:BP:V5:V6 and save this as reference for downstream operations
bim <- bim |>
  dplyr::mutate(
    snp_original = paste0(CHR, ":", BP, ":", V5, ":", V6)
  )

if (is_pgen) {
    pvar_standard <- data.table(
        `#CHROM` = bim$CHR,
        POS = bim$BP,
        ID = bim$snp_original,
        REF = bim$V5,
        ALT = bim$V6
    )
    fwrite(pvar_standard, paste0(opt$bfile, ".standard_snpid.pvar"), quote = FALSE, na = NA, sep = "\t")
} else {
    fwrite(bim, paste0(opt$bfile, ".standard_snpid.bim"), quote = FALSE, na = NA, sep = "\t", col.names = FALSE)
}



# Remove SNPs where SNP ID is duplicated since this can mess up liftOver and other operations downstream
if (sum(duplicated(bim$snp_original)) > 0 ) {
  message(sum(duplicated(bim$snp_original)), " duplicated snp ids in bim file. Removing them")
  bim <- bim[!(duplicated(bim[, .(snp_original)]) | duplicated(bim[, .(snp_original)], fromLast=T))]
}

# If necessary, lift to build 38
if(as.numeric(opt$grch)==37 && as.logical(opt$run_liftover)){
  bim_to_clean <- hg19ToHg38_liftover(bim)
} else {
  bim_to_clean <- bim
}

# Remove rows with duplicated SNP by CHR POS
# This get rid of multi-allelic variants and any other odd situations
bim_cleaned <- bim_to_clean[!(duplicated(bim_to_clean[, .(CHR, BP)]) | duplicated(bim_to_clean[, .(CHR, BP)], fromLast=T))]

# Save list of SNP ids to extract
extract_file <- paste0(opt$bfile, "_snps_to_extract.txt")
fwrite(list(bim_cleaned |> dplyr::pull(snp_original) |> unique()), extract_file, col.names=F, quote=F)

# Extract list of SNPs and output as bed/bim/fam or pgen/pvar/psam
if (is_pgen) {
  exit_status = system(paste0("plink2 --pgen ", opt$bfile, ".pgen --psam ", opt$bfile, ".psam --pvar ", opt$bfile, ".standard_snpid.pvar --extract ", extract_file, " --make-pgen --out ", opt$bfile, ".GRCh38.alpha_sorted_alleles"))
} else {
  exit_status = system(paste0("plink2 --bed ", opt$bfile, ".bed --fam ", opt$bfile, ".fam --bim ", opt$bfile, ".standard_snpid.bim --extract ", extract_file, " --make-bed --out ", opt$bfile, ".GRCh38.alpha_sorted_alleles"))
}
  
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

# Save
if (is_pgen) {
    pvar_alpha_sorted <- data.table(
      `#CHROM` = bim_alpha_sorted$CHR,
      POS = bim_alpha_sorted$BP,
      ID = bim_alpha_sorted$snp_original,
      REF = bim_alpha_sorted$V5,
      ALT = bim_alpha_sorted$V6
    )

    fwrite(
      pvar_alpha_sorted,
      paste0(opt$bfile, ".GRCh38.alpha_sorted_alleles.pvar"),
      quote = FALSE, na = NA, sep = "\t")
} else {
  fwrite(
    bim_alpha_sorted,
    paste0(opt$bfile, ".GRCh38.alpha_sorted_alleles.bim"),
    quote=F, na=NA, sep="\t", col.names = F)
}
