#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(readr))
suppressMessages(library(R.utils))
suppressMessages(library(data.table))
liftOver <- rtracklayer::liftOver
import.chain <- rtracklayer::import.chain
GRanges <- GenomicRanges::GRanges
IRanges <- IRanges::IRanges

### Munging function -------
#' Munge GWAS summary statistics
#'
#' This function processes a GWAS summary statistics file by renaming columns to standardized labels,
#' filtering for autosomal chromosomes (1-22), calculating MAF and sample size when missing,
#' and calculating the variance of effect size estimates if not provided by the user. 
#' Duplicated SNPs according to CHR:BP are removed to remove potential multi-allelic variants.
#' It supports both quantitative and case-control traits.
#'
#' @param sumstats.file A character string specifying the path to a GWAS summary statistics file, or a data frame containing the summary statistics.
#' @param snp.lab A character string specifying the name of the SNP column in the summary statistics dataset. Default is `"SNP"`.
#' @param chr.lab A character string specifying the name of the chromosome column. Default is `"CHR"`.
#' @param pos.lab A character string specifying the name of the base pair position column. Default is `"BP"`.
#' @param a1.lab A character string specifying the name of the effect allele column. Default is `"A1"`.
#' @param a0.lab A character string specifying the name of the non-effect allele column. Default is `"A2"`.
#' @param beta.lab A character string specifying the name of the effect size (beta) column. Default is `"BETA"`.
#' @param se.lab A character string specifying the name of the standard error column. Default is `"SE"`.
#' @param pval.lab A character string specifying the name of the p-value column. Default is `"P"`.
#' @param freq.lab A character string specifying the name of the allele frequency column. Default is `"FRQ"`.
#' @param n.lab A character string specifying the name of the sample size column. Default is `"N"`. If missing, it will be estimated from the allele frequency and standard error.
#' @param type A character string specifying the type of the trait. Options are `"quant"` for quantitative traits or `"cc"` for case-control traits. Default is `NULL`.
#' @param sdY A numeric value or character string specifying the standard deviation of the trait. For quantitative traits, if `sdY` is provided, it can be a single value or a file containing per-phenotype `sdY` values. If `NULL`, it will be estimated.
#' @param s A numeric value specifying the proportion of cases in case-control studies (required if `type` is `"cc"`). Default is `NULL`.
#'
#' @return A data frame containing the cleaned summary statistics, including standardized column names, and calculated or validated fields such as MAF, `varbeta`, and sample size (`N`).
#'
#' @details
#' - The function filters for autosomal chromosomes (1-22).
#' - If MAF is not provided, it is inferred from the allele frequency column (`freq`).
#' - For case-control studies, `s` (the proportion of cases) must be provided, while for quantitative traits, `sdY` can be estimated if not provided.
#' - The function can handle p-values reported in log10 scale, converting them to the original p-value.
#' - It removes any rows with missing values after all transformations are applied.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' munged_dataset <- dataset.munge_hor("gwas_sumstats.txt", type = "quant", sdY = 1.5)
#' }
#' 
#' @import dplyr
#' @import data.table
#' @export

dataset.munge_hor=function(sumstats.file
                           ,snp.lab="SNP"
                           ,chr.lab="CHR"
                           ,pos.lab="BP"
                           ,a1.lab="A1"
                           ,a0.lab="A2"
                           ,beta.lab="BETA"
                           ,se.lab="SE"
                           ,pval.lab="P"
                           ,freq.lab="FRQ"
                           ,n.lab="N"
                           ,type=NULL
                           ,sdY=NULL
                           ,s=NULL
                           ){
  
  # Load sumstat
  dataset = gwas
  if(is.character(sumstats.file)){
    dataset=read_delim(sumstats.file, data.table=F)
    dataset <- as.data.table(dataset)
  }
  
  if(!is.null(a1.lab) & a1.lab %in% names(dataset) & !is.null(a0.lab) & a0.lab %in% names(dataset) ){
    names(dataset)[match(c(a1.lab,a0.lab),names(dataset))]=c("A1","A2")
  }else{
    stop("a0.lab or a1.lab have not been defined or the column is missing")
  }
  if(!is.null(beta.lab) & beta.lab %in% names(dataset)){
    names(dataset)[names(dataset)==beta.lab]="BETA"
  }else{
    stop("beta.lab has not been defined or the column is missing")
  }
  
  if(!is.null(snp.lab) & snp.lab %in% names(dataset)){
    names(dataset)[names(dataset)==snp.lab]="snp_original"
  }else{
    stop("snp.lab has not been defined or the column is missing")
  }
  
  if(!is.null(se.lab) & se.lab %in% names(dataset)){
    names(dataset)[names(dataset)==se.lab]="SE"
  }else{
    stop("se.lab has not been defined or the column is missing")
  }
  
  if(!is.null(chr.lab) & chr.lab %in% names(dataset)){
    names(dataset)[names(dataset)==chr.lab]="CHR"
    dataset[ , CHR := as.numeric(dataset$CHR)] ### set as numeric, can't merge columns of different classes
    dataset <- dataset[CHR %in% c(1:22)]
  }else{
    stop("chr.lab has not been defined or the column is missing") ### Can be as well retrieved from LD reference bfiles??
  }
  
  if(!is.null(pos.lab) & pos.lab%in%names(dataset)){
    names(dataset)[names(dataset)==pos.lab]="BP"
    dataset[, BP := as.numeric(dataset$BP)] ### set as numeric, can't merge columns of different classes
  }else{
    stop("pos.lab has not been defined or the column is missing") ### Can be as well retrieved from LD reference bfiles??
  }
  
  #### Put either frequency of N mandatory
  if( !( freq.lab %in% names(dataset) | n.lab %in% names(dataset) )) {
    stop("Either effect allele frequency or N sample size needs to be provided!")
  }
  
  if(!is.null(freq.lab) & freq.lab %in% names(dataset)){
    names(dataset)[names(dataset)==freq.lab]="freq"
  } ### if effect allele frequency is not reported, it will calculated from the LD reference bfiles in the alignment step
  
  if("freq" %in% colnames(dataset)){
    dataset[, MAF := freq]
    dataset <- dataset[, MAF := ifelse(MAF<0.5, MAF, 1-MAF)]
  }
  
  if(!is.null(n.lab) & n.lab %in% names(dataset)){
    names(dataset)[names(dataset)==n.lab]="N"
  } else if( (is.null(n.lab) || is.na(n.lab) || n.lab=="NA" || !(n.lab %in% names(dataset))) && "MAF" %in% names(dataset)) {
    N_hat<-median(1/((2*dataset$MAF*(1-dataset$MAF))*dataset$SE^2),na.rm = T) 
    dataset$N=ceiling(N_hat)
  } ### if both effect allele, MAF and N are not reported, it will calculated from the LD reference bfiles in the alignment step
  
  if(!is.null(pval.lab) & pval.lab %in% names(dataset)){
    names(dataset)[names(dataset)==pval.lab]="P"
    #dataset$P <- as.numeric(dataset$P)
    ### Remove NA p-values
    dataset <- dataset[!is.na(P)]
    ### Check if p-value column provided is log10 transformed. If yes, compute original p-value
    if (!all(dataset$P >= 0 & dataset$P <= 1)) {
      dataset <- dataset[, P:=10^(-P)]
    }
  } else {
    dataset[, P := pchisq((BETA/SE)^2,df=1,lower=F)]
  }
  # Add variance of beta  
  dataset[, varbeta := SE^2]
  
  # Add type and sdY/s
  dataset[, type:= type]
  if(type=="cc" & !(is.null(s)) && !is.na(s)){ # && prevents to return "logical(0)" when s is null
    dataset[, s := s]
    colname_for_type <- "s"
  } else if(type=="cc" & (is.null(s) || is.na(s))){
    #### Is this correct?? Is "s" strictly necessary for cc traits??
    stop("Please provide s, the proportion of samples who are cases")
  }
  
  if (type == "quant") {
    colname_for_type <- "sdY"
    # If multiple values (one per key if the GWAS is mol_QTL) are provided in a file
    if (is.character(sdY) && file.exists(sdY)) {
      sdY_list <- read_delim(sdY, data.table = F)
      sdY_list <- as.data.table(sdY_list)
      
      if(all(c("sdY", "phenotype_id") %in% names(sdY_list))){
        if(is.numeric(sdY_list$sdY) & all(dataset$phenotype_id %in% sdY_list$phenotype_id)){
          dataset <- merge(dataset, sdY_list, by = "phenotype_id", all.x=T)
        }
      } else {
        stop("If you're providing sdY through a file check 1) that the file exists, 2) that it has a \"sdY\" and \"phenotype_id\" column and 3) that a sdY is provided for ALL traits in the GWAS sum stat")
      }
      
      # If a single value is provided
    } else if (!is.null(sdY) && !is.na(sdY) && sdY != "NA") {
      dataset[ , sdY := sdY]
    } else if ("MAF" %in% names(dataset)) {
      # If sdY is not provided, calculate it
      dataset[, sdY := coloc:::sdY.est(varbeta, MAF, N), by = phenotype_id]
    }
  }
  
  # Make a temporary SNPID column as CHR:BP as remove duplicated SNPs
  dataset[, TMP_SNPID := paste0("chr", CHR, ":", BP)]
  n_dups_removed <- sum(duplicated(dataset$TMP_SNPID) | duplicated(dataset$TMP_SNPID, fromLast = TRUE))
  dataset <- dataset[!duplicated(TMP_SNPID) | !duplicated(TMP_SNPID, fromLast = TRUE)]
  if (n_dups_removed > 0) {
    message("Removed ", n_dups_removed, " duplicated SNPs based on CHR:BP")
  }
  
  # Select only necessary columns
  columns <- c("phenotype_id","snp_original","CHR","BP","A1","A2","BETA","varbeta","SE","P","type", intersect(c("N","freq","MAF","s", "sdY"), names(dataset)))
  cols_to_remove <- setdiff(names(dataset), columns)
  dataset[, (cols_to_remove) := NULL]
  
  if(type == "quant" && "s" %in% names(dataset)){
    dataset$s <- NULL
  } else if(type == "cc" && "sdY" %in% names(dataset)) {
    dataset$sdY <- NULL
  }
  
  # Remove all rows with NAs 
  dataset <- dataset[complete.cases(dataset)] # Remove rows with any NA values.
  setorder(dataset, CHR, BP) # Sort by CHR then BP
  
  return(dataset)
}

### liftover function ---------
hg19ToHg38_liftover <- function(
    dataset_munged,
    default_chain_file = "hg19ToHg38.over.chain"
) {
  
  chain_file <- default_chain_file
  
  if (!file.exists(chain_file)) {
    message("Downloading chain file...")
    download.file("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                  destfile = paste0(chain_file, ".gz"))
    R.utils::gunzip(paste0(chain_file, ".gz"))
  }
  
  ch <- import.chain(chain_file)
  
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


### dataset.align function --------
dataset.align <- function(dataset, bfile) {
  
  # Load munged dataset sumstat in .rds format (if necessary)
  if(is.character(dataset)){
    dataset <- read_delim(dataset, data.table=F)
    dataset <- as.data.table(dataset)
  }
  
  setnames(dataset, c("A1", "A2"), c("A1_old", "A2_old")) # Rename A1 and A2
  
  dataset[, `:=`(
    A1 = pmin(A1_old, A2_old),
    A2 = pmax(A1_old, A2_old)
  )]
  
  dataset[, SNP := sprintf("chr%s:%s:%s:%s", CHR, BP, A1, A2)]
  dataset[, b := fcase(
    A1_old == A2 & A2_old == A1, -BETA,
    A1_old == A1 & A2_old == A2, BETA 
  )]
  
  ### Check if freq (and MAF) are present in the dataframe - if not, compute them from LD reference
  
  if("freq" %in% colnames(dataset)){
    setnames(dataset, "freq", "freq_old") # Rename 'freq' to 'freq_old'
    dataset[, freq := ifelse(A1_old == A2 & A2_old == A1, (1 - freq_old), freq_old)] # Calculate new 'freq'
    dataset[, freq_old := NULL] # Remove the old column
  } else {
    message("Allele frequency not found in summary stat - Computing allele frequency from LD reference panel")
    ## Compute allele frequency from LD reference panel provided    
    random.number <- stringi::stri_rand_strings(n=1, length=20, pattern = "[A-Za-z0-9]")
    exit_status = system(paste0("plink2 --bfile ", bfile, " --freq --make-just-bim --out ", random.number))

    # Raise an error if the external command fails
    if (exit_status != 0) {
      cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
      quit(status = 1, save = "no")
    }

    # Load-in frequency and position info from plink files
    # freqs <- as.data.table(cbind(
    #   read_delim(paste0(random.number,".bim"), col_names = FALSE) %>% dplyr::select(X1, X4) %>% dplyr::rename(CHR=X1, BP=X4),
    #   read_delim(paste0(random.number,".afreq")) %>% dplyr::select(REF, ALT, ALT_FREQS)
    # ))
    
    # Read .bim and select/rename columns
    bim <- read_delim(paste0(random.number, ".bim"), col_names = FALSE)
    bim <- bim[, c(1, 4)]
    names(bim) <- c("CHR", "BP")
    
    # Read .afreq and select needed columns
    afreq <- read_delim(paste0(random.number, ".afreq"))
    afreq <- afreq[, c("REF", "ALT", "ALT_FREQS")]
    
    # Combine and convert to data.table
    freqs <- as.data.table(cbind(bim, afreq))
    
    # Compute frequency for effect allele, create SNP column for merging
    freqs[, `:=`(
      A1 = pmin(REF, ALT),
      A2 = pmax(REF, ALT)
    )]
    
    freqs[, freq := ifelse(ALT == A1 & REF == A2, ALT_FREQS, 1 - ALT_FREQS)]
    
    freqs[, `:=`(
      MAF = ifelse(freq < 0.5, freq, 1 - freq),
      SNP = paste0("chr", CHR, ":", BP, ":", A1, ":", A2)
    )]
    
    freqs <- freqs[, .(SNP, freq, MAF)]
    
    # Add freq and MAF info to dataset 
    dataset <- merge(dataset, freqs, by = "SNP")
    
    system(paste0("rm ",random.number,"*"))
  }
  
  # Remove the old columns
  dataset[, `:=`(A1_old = NULL, A2_old = NULL)]
  
  #### You can finally calculate N!
  if(!("N" %in% colnames(dataset))){
    N_hat <- median(1/((2*dataset$MAF*(1-dataset$MAF))*dataset$SE^2),na.rm = T) 
    dataset[ , N := ceiling(N_hat)]
  }
  
  #### You can finally calculate sdY! If sdY is not provided, calculate it - ONLY IF type is not cc or sdY has not been already provided or calculated
  if( !(unique(dataset$type)=="cc" | "sdY" %in% names(dataset)) ){
    dataset[, sdY := coloc:::sdY.est(varbeta, MAF, N), by = phenotype_id]
  }
  
  # Select columns in the specified order
  cols <- c("snp_original", "SNP", "CHR", "BP", "A1", "A2", "freq", "b", "varbeta", "SE", "P", "MAF", "N", "type", intersect(c("s", "sdY"), names(dataset)), "phenotype_id")
  cols_to_remove <- setdiff(names(dataset), cols)
  dataset[, (cols_to_remove) := NULL]
  
  # Rename columns
  setnames(dataset, c("SE", "P"), c("se", "p"))
  setkey(dataset, SNP, phenotype_id)
  
  # Remove all rows with duplicated values in the SNP column, by phenotype
  # This is needed to ensure that reordering of alleles does not create duplicated SNPs
  dataset <- dataset[, .SD[!duplicated(SNP) & !duplicated(SNP, fromLast =TRUE)], by = phenotype_id]
  
  return(dataset)
}


### locus.breaker function --------
locus.breaker <- function(
    res,
    p.sig = 5e-08,
    p.limit = 1e-05,
    hole.size = 250000,
    p.label = "P",
    chr.label = "CHR",
    pos.label = "BP"){
  
  #res <- as.data.frame(res)
  setorderv(res, c(chr.label, pos.label))
  res = res[get(p.label) < p.limit]
  trait.res = c()
  
  for(j in unique(res[, get(chr.label)])) {
    res.chr = res[get(chr.label) == j]
    if (nrow(res.chr) > 1) {
      pos_values <- res.chr[, get(pos.label)] # Get the pos.label column values
      holes <- pos_values[-1] - pos_values[-length(pos_values)]
      gaps = which(holes > hole.size)
      if (length(gaps) > 0) {
        for (k in 1:(length(gaps) + 1)) {
          if (k == 1) {
            res.loc = res.chr[1:(gaps[k]), ]
          }
          else if (k == (length(gaps) + 1)) {
            res.loc = res.chr[(gaps[k - 1] + 1):nrow(res.chr), 
            ]
          } else {
            res.loc = res.chr[(gaps[k - 1] + 1):(gaps[k]), 
            ]
          }
          if (min(res.loc[, get(p.label)]) < p.sig) {
            start.pos = min(res.loc[, get(pos.label)], na.rm = T)
            end.pos = max(res.loc[, get(pos.label)], na.rm = T)
            chr = j
            best.snp = res.loc[which.min(res.loc[, get(p.label)]), 
            ]
            line.res = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (min(res.loc[, get(p.label)]) < p.sig) {
          start.pos = min(res.loc[, get(pos.label)], na.rm = T)
          end.pos = max(res.loc[, get(pos.label)], na.rm = T)
          chr = j
          best.snp = res.loc[which.min(res.loc[, get(p.label)]), 
          ]
          line.res = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (min(res.loc[, get(p.label)]) < p.sig) {
        start.pos = min(res.loc[, get(pos.label)], na.rm = T)
        end.pos = max(res.loc[, get(pos.label)], na.rm = T)
        chr = j
        best.snp = res.loc[which.min(res.loc[,get(p.label)]), 
        ]
        line.res = c(chr, start.pos, end.pos, unlist(best.snp))
        trait.res = rbind(trait.res, line.res)
      }
    }
  }
  if(!is.null(trait.res)){
    trait.res = as.data.table(trait.res, stringsAsFactors = FALSE)
    trait.res[, (chr.label) := NULL] # Remove the column specified by chr.label
    setnames(trait.res, 1:3, c("chr", "start", "end")) # Rename the first 3 columns
    setattr(trait.res, "row.names", NULL) # Remove row names
  }
  return(trait.res)
}

### Accessory functions --------
read_first_line <- function(file_path, sep="\t") {
  if (endsWith(file_path, ".gz")) {
    con <- gzfile(file_path, "rt")
  } else {
    con <- file(file_path, "rt")
  }
  
  on.exit(close(con)) # Ensure connection is closed
  
  first_line <- readLines(con, n = 1)
  if (length(first_line) == 0) {
    return(NULL)
  }
  elements <- strsplit(first_line, sep)[[1]]
  
  return(elements)
}

find_positions <- function(A, B) {
  positions_A_in_B <- match(A, B)
  positions_A_in_B <- positions_A_in_B[!is.na(positions_A_in_B)]
  
  return(positions_A_in_B)
}

# Function to round scientific notation to a specific number of decimal places
# If we have less than 17 decimals we return the number as is
round_sci <- function(x, non_sci_digits = 17, sci_digits=15) {
  x_str <- format(x, scientific = FALSE)
  decimal_part <- sub("^[^.]*\\.", "", x_str)
  zeros <- regmatches(decimal_part, regexpr("^0*", decimal_part))

  if (nchar(zeros) < 4) {
    formatted_value <- format(round(x,non_sci_digits), nsmall=non_sci_digits, scientific = FALSE)
  } else {
    formatted_value <- format(x, digits=sci_digits, scientific = TRUE)
  }
  return(formatted_value)
}

# Get arguments --------
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
  make_option("--sdY", default=NULL, help="For a quantitative trait (type==quant), the population standard deviation of the trait. For quantitative traits, it can be a single value or a file containing per-phenotype `sdY` values. If not given, it will be estimated from beta and MAF"),
  make_option("--s", default=NULL, help="For a case control study (type==cc), the proportion of samples in dataset 1 that are cases"),
  make_option("--bfile", default=NULL, help="Path and prefix name of custom/default LD bfiles (PLINK format .bed .bim .fam) - to compute effect allele frequency if missing"),
  make_option("--grch", default=NULL, help="Genome reference build of GWAS sum stats"),
  make_option("--run_liftover", default=FALSE, type="logical", help="Set true to run liftover"),
  make_option("--maf", default=1e-04, help="MAF filter", metavar="character"),
  make_option("--p_thresh1", default=5e-08, help="Significant p-value threshold for top hits"),
  make_option("--p_thresh2", default=1e-05, help="P-value threshold for loci borders"),  
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--study_id", default=NULL, help="Id of the study"),
  make_option("--threads", default=1, help="N threads for processing")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
#source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once_vroom.R"))

## Throw error message - GWAS summary statistics file MUST be provided!
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Please specify the path and file name of your GWAS summary statistics in --path option", call.=FALSE)
}


###### NB: Move all checks (files correctly provided, reasonable values etc.) in a specific function!!


##################################
# Load in and munge GWAS sum stats
##################################

message("Reading data from ", opt$input)
input_colnames <- read_first_line(opt$input)
integer_columns <- c(opt$pos_lab, opt$n_lab)
double_columns <- c(opt$freq_lab, opt$effect_lab, opt$se_lab, opt$pvalue_lab, opt$sdY)

idx_int_columns<- find_positions(integer_columns, input_colnames)
idx_dbl_columns <- find_positions(double_columns, input_colnames)

dtypes <- rep("c", length(input_colnames))
dtypes[idx_int_columns] <- 'i'
dtypes[idx_dbl_columns] <- 'd'
dtypes <- paste(dtypes, collapse="")

gwas <- read_delim(opt$input, na = c("", "NA"), num_threads = opt$threads, col_types = dtypes, lazy=TRUE) # treat both "NA" as character and empty strings ("") as NA
gwas <- as.data.table(gwas)
#gwas <- fread(opt$input, na.strings = c("", "NA"), tmpdir=getwd()) # treat both "NA" as character and empty strings ("") as NA

# If the trait column is NOT provided, add the same one for the whole sum stat
if(isFALSE(opt$is_molQTL)){
  gwas[, phenotype_id := "full"]
  opt$key="phenotype_id"
# If the trait column is provided, simply rename it
} else if (isTRUE(opt$is_molQTL)) {
  setnames(gwas, opt$key, "phenotype_id")
}

gc()

# Format GWAS summary statistics
message("Munging")
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
rm(gwas)
gc()

# If necessary, lift to build 38
if(as.numeric(opt$grch)==37 & isTRUE(opt$run_liftover)){
  message("Performing liftOver to GRCh38")
  # To ensure coordinate mapping works correctly, we temporarly set snp_original to CHR:BP:A1:A2
  # We store the snp_original value in a temporary column
  dataset_munged[, TMP_SNPID := snp_original] 
  dataset_munged[, snp_original := paste0("chr", CHR, ":", BP, ":", A1, ":", A2)] 

  dataset_munged <- hg19ToHg38_liftover(dataset_munged)

  # Restore the original snp_original column
  dataset_munged[, snp_original := TMP_SNPID]
  dataset_munged[, TMP_SNPID := NULL] # Remove the temporary column
  gc()
}

# Align GWAS sum stats
message("Align alleles to BIM file")
dataset_munged <- dataset.align(dataset_munged, bfile=opt$bfile)

# Now that we have clean, aligned data, do a final check for duplicated CHR BP and remove all duplicates
# Report the number of duplicates removed
n_vars_before <- nrow(dataset_munged)
dataset_munged <- dataset_munged[, .SD[!(duplicated(.SD[, .(CHR, BP)]) | duplicated(.SD[, .(CHR, BP)], fromLast = TRUE))], by = phenotype_id]
n_vars_removed <- n_vars_before - nrow(dataset_munged)
if (n_vars_removed > 0) {
  message("Removed ", n_vars_removed, " during final check of duplicated CHR:BP variants")
}

# Perform MAF filter here (doesn't make sense to do this later!)
dataset_munged <- dataset_munged[MAF > opt$maf]
setorder(dataset_munged, CHR) # Sort by CHR
gc()

################
# Locus breaker
################

message("Run LOCUS BREAKER")
message("Columns in input data: ", paste(names(dataset_munged), collapse=";"))

loci_list <- dataset_munged[, {
  if (sum(.SD$p < opt$p_thresh1) > 0) {
    locus.breaker(
      .SD,
      p.sig = opt$p_thresh1,
      p.limit = opt$p_thresh2,
      hole.size = opt$hole,
      p.label = "p",
      chr.label = "CHR",
      pos.label = "BP"
    )
  }
}, by = phenotype_id]

# Slightly enlarge locus by 200kb!
if(nrow(loci_list) > 0){
  loci_list[, start := as.numeric(start) - 100000]
  loci_list[start < 0, start := 0] # Replace negative numbers with 0
  loci_list[, end := as.numeric(end) + 100000]
  
  ### Add study ID to the loci table. Save
  loci_list[, study_id := opt$study_id]
  
  ### Check if locus spans the HLA locus chr6:28,510,120-33,480,577
  ### https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
  # This code checks for four conditions to determine if a locus partially or completely spans the HLA region:
  # 1) Locus starts before (or at the exact beginning of) HLA and ends after (or at the exact end of) HLA.
  # 2) Locus starts and ends within the HLA region.
  # 3) Locus starts before the end of HLA and ends after the end of HLA.
  # 4) Locus starts before the beginning of HLA and ends after the beginning of HLA.
  hla_start=28510120
  hla_end=33480577
  hla_coord <- seq(hla_start,hla_end)
  
  loci_list[, is_in_hla := chr == 6 & 
              ((start <= hla_end & end >= hla_start) | 
                 (start >= hla_start & end <= hla_end) | 
                 (start <= hla_end & end >= hla_end) | 
                 (start <= hla_start & end >= hla_start))]
  loci_list[, locus_size := end - start]

  cat(paste0("\n", nrow(loci_list), " significant loci identified for ", opt$study_id, "\n"))
  cat(paste0("Loci size range: ", min(loci_list$locus_size), " - ", max(loci_list$locus_size), "\n"))
  
  # Reorder columns in the loci_list table
  columns_order <- c("chr", "start", "end", "locus_size", "phenotype_id", "snp_original", "SNP", "BP", "A1", "A2", "freq", "b", "varbeta", "se", "p", "MAF", "N", "type", "s", "sdY", "study_id", "is_in_hla")
  columns_order <- intersect(columns_order, names(loci_list))
  message("Final column order: ", paste(columns_order, collapse=","))
  setcolorder(loci_list, columns_order)
  
  write_tsv(loci_list, paste0(opt$study_id, "_loci.tsv"), num_threads = opt$threads, quote="none", na="NA") ### full table to publish
  write_tsv(loci_list[is_in_hla=="FALSE"], paste0(opt$study_id, "_loci_NO_HLA.tsv"), num_threads = opt$threads, quote="none", na="NA") ### table to feed to fine-mapping - for the moment, HLA region ALWAYS removed
  
} else {
  message("No loci detected - skip")
}
gc()

message("Save processed dataset")
### Order by phenotype_id, which will be used by tabix to index
# Reorder columns in the output table
columns_order <- c("phenotype_id", "snp_original", "SNP", "CHR", "BP", "A1", "A2", "freq", "b", "se", "p", "N", "type", "sdY", "s")
columns_order <- intersect(columns_order, names(dataset_munged))
setcolorder(dataset_munged, columns_order)
setorder(dataset_munged, phenotype_id, BP)

# round float to a fixed precision
dataset_munged[, freq := round(freq, 7)]
dataset_munged[, p := sapply(p, function(x) round_sci(x))]

### Save removing info you don't need - varbeta is later calculated on bC_se
write_tsv(dataset_munged[, `:=`(MAF = NULL, varbeta = NULL)], paste0(opt$study_id, "_dataset_aligned.tsv.gz"), num_threads = opt$threads, quote="none", na="NA")

### BGZIP compress and index the file
from <- paste0(opt$study_id, "_dataset_aligned.tsv.gz")
to <- paste0(opt$study_id, "_dataset_aligned_indexed.tsv.gz")
zipped <- Rsamtools::bgzip(from, to, overwrite = T)
idx <- Rsamtools::indexTabix(zipped,
                  skip=as.integer(1),
                  seq=which(colnames(dataset_munged)=="phenotype_id"),
                  start=which(colnames(dataset_munged)=="BP"),
                  end=which(colnames(dataset_munged)=="BP")
                  )

message("-- COMPLETED SUCCESSFULLY --")