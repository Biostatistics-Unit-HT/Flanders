suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--input", default=NULL, help="Path and file name of GWAS summary statistics"),
  make_option("--is_molQTL", default=NULL, type="logical", help="Whether the summary statistics provided are for eQTLs, ATAC, ChipSeq etc. data"),
  make_option("--key", default=NULL, help="If GWAS is from a molQTL, name of column reporting the trait/phenoptype"),
  make_option("--chr_lab", default="CHROM", help="Name of chromosome column in GWAS summary statistics"),
  make_option("--pos_lab", default="GENPOS", help="Name of genetic position column in GWAS summary statistics"),
  make_option("--rsid_lab", default="ID", help="Name of rsid column"),  
  make_option("--a1_lab", default="ALLELE1", help="Name of effect allele column"),  
  make_option("--a0_lab", default="ALLELE0", help="Name of NON effect allele column"),  
  make_option("--freq_lab", default="A1FREQ", help="Name of effect allele frequency column"),  
  make_option("--n_lab", default="N", help="Name of sample size column"),  
  make_option("--effect_lab", default="BETA", help="Name of effect size column"),  
  make_option("--se_lab", default="SE", help="Name of standard error of effect column"), 
  make_option("--pvalue_lab", default="P", help="Name of p-value of effect column"),
  make_option("--type", default=NULL, help="Type of phenotype analysed - either 'quant' or 'cc' to denote quantitative or case-control"),
  make_option("--sdY", default=NULL, help="For a quantitative trait (type==quant), the population standard deviation of the trait. If not given, it will be estimated beta and MAF"),
  make_option("--s", default=NULL, help="For a case control study (type==cc), the proportion of samples in dataset 1 that are cases"),
  make_option("--bfile", default=NULL, help="Path and prefix name of custom/default LD bfiles (PLINK format .bed .bim .fam) - to compute effect allele frequency if missing"),
  make_option("--grch", default=NULL, help="Genome reference build of GWAS sum stats"),
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

## Throw error message - GWAS summary statistics file MUST be provided!
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Please specify the path and file name of your GWAS summary statistics in --path option", call.=FALSE)
}


###### NB: Move all checks (files correctly provided, reasonable values etc.) in a specific function!!

##################################
# Check columns of GWAS
##################################

# Read only the header (no rows) of the file
gwas_header = fread(opt$input, na.strings = c("", "NA"), tmpdir=getwd(), nrow = 0)
required_columns <- c(opt$rsid_lab, opt$chr_lab, opt$pos_lab, opt$a1_lab, opt$a0_lab, opt$effect_lab, opt$se_lab)
missing_columns <- setdiff(required_columns, colnames(gwas_header))
if (length(missing_columns) > 0) {
  stop(paste0("Error: The following required columns are missing in the GWAS file ", opt$input, ": ", paste(missing_columns, collapse = ", ")))
}

##################################
# Load in and munge GWAS sum stats
##################################

gwas <- fread(opt$input, na.strings = c("", "NA"), tmpdir=getwd()) # treat both "NA" as character and empty strings ("") as NA

# If the trait column is NOT provided, add the same one for the whole sum stat
if(isFALSE(opt$is_molQTL)){
  gwas <- gwas %>% mutate(phenotype_id="full")
  opt$key="phenotype_id"
# If the trait column is provided, simply rename it
} else if (isTRUE(opt$is_molQTL)) {
  gwas <- gwas %>% rename(phenotype_id=opt$key)
}

# Format GWAS summary statistics
dataset_munged <- dataset.munge_hor(
  gwas
  ,snp.lab = opt$rsid_lab
  ,chr.lab = opt$chr_lab
  ,pos.lab = opt$pos_lab
  ,a1.lab = opt$a1_lab
  ,a0.lab = opt$a0_lab
  ,beta.lab = opt$effect_lab
  ,se.lab = opt$se_lab
  ,freq.lab = opt$freq_lab
  ,pval.lab = opt$pvalue_lab
  ,n.lab = opt$n_lab
  ,type = opt$type
  ,sdY = opt$sdY
  ,s = opt$s
)

# If necessary, lift to build 38
if(as.numeric(opt$grch)==37){
  dataset_munged <- hg19ToHg38_liftover(dataset_munged)
}

# Align GWAS sum stats
dataset_aligned <- dataset.align(dataset_munged, bfile=opt$bfile)

# Perform MAF filter here (doesn't make sense to do this later!)
dataset_aligned <- dataset_aligned %>%
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
  fwrite(loci_list, paste0(opt$study_id, "_loci.tsv"), sep="\t", quote=F, na=NA) ### full table to publish
  fwrite(loci_list %>% dplyr::filter(is_in_hla=="FALSE"), paste0(opt$study_id, "_loci_NO_HLA.tsv"), sep="\t", quote=F, na=NA) ### table to feed to fine-mapping - for the moment, HLA region ALWAYS removed
  
}

### Order by phenotype_id, which will be used by tabix to index
dataset_aligned <- as.data.table(dataset_aligned)
dataset_aligned$BP <- as.integer(dataset_aligned$BP)
setorder(dataset_aligned, phenotype_id, BP)

### Save removing info you don't need - varbeta is later calculate on bC_se
fwrite(dataset_aligned %>% select(-MAF, -varbeta), paste0(opt$study_id, "_dataset_aligned.tsv.gz"), sep="\t", quote=F, na=NA)

#from <- dataset_aligned ### way slower setting from to the R object rather than file
from <- paste0(opt$study_id, "_dataset_aligned.tsv.gz")
to <- paste0(opt$study_id, "_dataset_aligned_indexed.gz")

### Convert .gz in bgzip
zipped <- Rsamtools::bgzip(from, to)

### Index the File with tabix
idx <- Rsamtools::indexTabix(zipped,
                  skip=as.integer(1),
                  seq=which(colnames(dataset_aligned)=="phenotype_id"),
                  start=which(colnames(dataset_aligned)=="BP"),
                  end=which(colnames(dataset_aligned)=="BP")
                  )
