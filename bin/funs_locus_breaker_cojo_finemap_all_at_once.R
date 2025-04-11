# Load packages
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(corrplot))
suppressMessages(library(coloc))
suppressMessages(library(susieR))
suppressMessages(library(bigsnpr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(stringi))
suppressMessages(library(stringr))
suppressMessages(library(patchwork))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(igraph))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(plyr))
suppressMessages(library(Gviz))
suppressMessages(library(Matrix))
#suppressMessages(library(Rfast))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
#suppressMessages(library(EnsDb.Hsapiens.v86))
#suppressMessages(library(ggtext)) ### TO ADD TO CONDA ENV!
suppressMessages(library(dplyr))


#' Munge GWAS summary statistics
#'
#' This function processes a GWAS summary statistics file by renaming columns to standardized labels,
#' filtering for autosomal chromosomes (1-22), calculating MAF and sample size when missing,
#' and calculating the variance of effect size estimates if not provided by the user. It supports both quantitative and case-control traits.
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
  if(is.character(sumstats.file)){
    dataset=fread(sumstats.file, data.table=F)
  }else{
    dataset=as.data.frame(sumstats.file)
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
    dataset$CHR <- as.numeric(dataset$CHR) ### set as numeric, can't merge columns of different classes
    dataset <- dataset %>% dplyr::filter(CHR %in% c(1:22))
  }else{
    stop("chr.lab has not been defined or the column is missing") ### Can be as well retrieved from LD reference bfiles??
  }
  
  if(!is.null(pos.lab) & pos.lab%in%names(dataset)){
    names(dataset)[names(dataset)==pos.lab]="BP"
    dataset$BP <- as.numeric(dataset$BP) ### set as numeric, can't merge columns of different classes
  }else{
    stop("pos.lab has not been defined or the column is missing") ### Can be as well retrieved from LD reference bfiles??
  }
  
  if(!is.null(freq.lab) & freq.lab %in% names(dataset)){
    names(dataset)[names(dataset)==freq.lab]="freq"
  } ### if effect allele frequency is not reported, it will calculated from the LD reference bfiles in the alignment step

  if("freq" %in% colnames(dataset)){
    dataset$MAF=dataset$freq
    dataset <- dataset %>% mutate(MAF=ifelse(MAF<0.5, MAF, 1-MAF))
  }
  
  if (!is.null(n.lab) && !is.na(n.lab) && n.lab %in% names(dataset) && !all(is.na(dataset[[n.lab]]))) {
    names(dataset)[names(dataset)==n.lab]="N"
  } else if( "MAF" %in% names(dataset)) {
    N_hat <- median(1/((2*dataset$MAF*(1-dataset$MAF))*dataset$SE^2),na.rm = T) 
    dataset$N=ceiling(N_hat)
  } ### if both effect allele, MAF and N are not reported, it will calculated from the LD reference bfiles in the alignment step
  
  if(!is.null(pval.lab) & pval.lab %in% names(dataset)){
    names(dataset)[names(dataset)==pval.lab]="P"
    dataset$P <- as.numeric(dataset$P)
### Remove NA p-values
    dataset <- dataset %>% dplyr::filter(!is.na(P))
    ### Check if p-value column provided is log10 transformed. If yes, compute original p-value
    if (!all(dataset$P >= 0 & dataset$P <= 1)) {
      dataset <- dataset %>% mutate(P=10^(-P))
    }
  } else {
    dataset$P=pchisq((dataset$BETA/dataset$SE)^2,df=1,lower=F)
  }
  # Add variance of beta  
  dataset$varbeta=dataset$SE^2
  
  # Add type and sdY/s
  dataset$type <- type
  if(type=="cc" & !(is.null(s)) && !is.na(s)){ # && prevents to return "logical(0)" when s is null
    dataset$s=s
  } else if(type=="cc" & (is.null(s) || is.na(s))){
    #### Is this correct?? Is "s" strictly necessary for cc traits??
    stop("Please provide s, the proportion of samples who are cases")
  }
  
  if (type == "quant") {
# If multiple values (one per key if the GWAS is mol_QTL) are provided in a file
    if (is.character(sdY) && file.exists(sdY)) {
      sdY_list <- fread(sdY, data.table = F)

      if(all(c("sdY", "phenotype_id") %in% colnames(sdY_list))){
        if(is.numeric(sdY_list$sdY) & all(dataset$phenotype_id %in% sdY_list$phenotype_id)){
          dataset <- dataset %>% left_join(sdY_list, by = "phenotype_id")
        }
      } else {
        stop("If you're providing sdY through a file check 1) that the file exists, 2) that it has a \"sdY\" and \"phenotype_id\" column and 3) that a sdY is provided for ALL traits in the GWAS sum stat")
      }
      
# If a single value is provided
    } else if (!is.null(sdY) && !is.na(sdY) && sdY != "NA") {
      dataset$sdY <- sdY
    } else {
# If sdY is not provided, calculate it
      dataset <- dataset %>%
        dplyr::group_by(phenotype_id) %>%
        dplyr::mutate(sdY = coloc:::sdY.est(varbeta, MAF, N)) %>%
        as.data.frame()
    }
  }

# Select only necessary columns
  dataset <- dataset %>%
    dplyr::select(any_of(c("phenotype_id","snp_original","CHR","BP","A1","A2","freq","BETA","varbeta","SE","P","MAF","N","type", "s", "sdY")))
    
  if(type == "quant" && "s" %in% names(dataset)){
    dataset <- dataset %>% dplyr::select(-"s")
  } else if(type == "cc" && "sdY" %in% names(dataset)) {
    dataset <- dataset %>% dplyr::select(-"sdY")
  }

# Remove all rows with NAs 
  dataset <- dataset %>% na.omit() %>% dplyr::arrange(CHR, BP)
  dataset
}



### from hg19 To Hg38 liftover function
hg19ToHg38_liftover <- function(
    dataset_munged,
    default_ch_path="/processing_data/reference_datasets/liftOver/2023.1/Homo_sapiens/hg19/"
){
  
  if(!(file.exists(paste0(default_ch_path, "hg19ToHg38.over.chain")))){
    ### Download chain file and unzip it
    exit_status = system("wget -O - http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz | gunzip -c > hg19ToHg38.over.chain")

    # Raise an error if the external command fails
    if (exit_status != 0) {
      cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
      quit(status = 1, save = "no")
    }
    
    # Import from hg19 To Hg38 chain and then delete it
    ch <- import.chain("hg19ToHg38.over.chain")
    system("rm hg19ToHg38.over.chain")
  } else {
    # Import from hg19 To Hg38 chain already present on hpc
    ch <- import.chain(paste0(default_ch_path, "hg19ToHg38.over.chain"))
  }
  
  # Convert to GRanges object (assuming GENPOS is 1-based)
  dataset_ranges <- makeGRangesFromDataFrame(
    dataset_munged %>% dplyr::mutate(start=BP) %>% dplyr::distinct(snp_original, CHR, start, BP),
    start.field="start",
    end.field="BP",
    seqnames.field="CHR",
    keep.extra.columns=TRUE
  )
  dataset_ranges@seqnames@values <- paste0("chr", dataset_ranges@seqnames@values)  
  
  # Perform liftover
  dataset_ranges38 <- liftOver(dataset_ranges, ch)
  dataset_ranges38_df <- as.data.frame(unlist(dataset_ranges38)) %>%
    dplyr::rename(BP=end) %>%
#    dplyr::rename(BP=start) %>%
    dplyr::select(BP, snp_original)
  
  # Add lifted positions to GWAS sum stats
  dataset_lifted <- dataset_munged %>%
    dplyr::select(-BP) %>%
    inner_join(dataset_ranges38_df, by="snp_original")
  
  return(dataset_lifted)
}



### dataset.align ###
dataset.align <- function(dataset,
                          bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/p01_output/ukbb_all_30000_random_unrelated_white_british"
                          ){
  
  # Load munged dataset sumstat in .rds format (if necessary)
  if(is.character(dataset)){
    dataset_munged <- fread(dataset, data.table=F)
  }else{
    dataset_munged <- as.data.frame(dataset)
  }
  
  
  df_sorted <- dataset_munged %>%
    dplyr::rename(A1_old=A1, A2_old=A2) %>%
    dplyr::mutate(
      A1 = pmin(A1_old, A2_old), # Sort A1 and A2 alphabetically
      A2 = pmax(A1_old, A2_old),
      b = ifelse(A1_old == A2 & A2_old == A1, -BETA, BETA), # Change sign of BETA if order changed
      SNP = paste0("chr", CHR, ":", BP, ":", A1, ":", A2)
    )
  
### Check if freq (and MAF) are present in the dataframe - if not, compute them from defaul/custom LD reference

  if("freq" %in% colnames(df_sorted)){
    df_final <- df_sorted %>%
      dplyr::rename(freq_old=freq) %>%
      dplyr::mutate(freq = ifelse(A1_old == A2 & A2_old == A1, (1-freq_old), freq_old)) # Reverse frequency if order changed!
  
  } else {

## Compute allele frequency from LD reference panel provided    
    random.number <- stri_rand_strings(n=1, length=20, pattern = "[A-Za-z0-9]")
    exist_status = system(paste0("plink2 --bfile ", bfile, " --freq --make-bed --out ", random.number))

    # Raise an error if the external command fails
    if (exit_status != 0) {
      cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
      quit(status = 1, save = "no")
    }

      # Load-in frequency and position info from plink files
    freqs <- cbind(
      fread(paste0(random.number,".bim")) %>% dplyr::select(V1, V4) %>% dplyr::rename(CHR=V1, BP=V4),
      fread(paste0(random.number,".afreq")) %>% dplyr::select(REF, ALT, ALT_FREQS)
    )
        
# Compute frequency for effect allele, create SNP column for merging
    freqs <- freqs %>%
      dplyr::mutate(
        A1 = pmin(REF, ALT), # Sort A1 and A2 alphabetically
        A2 = pmax(REF, ALT),
        freq = ifelse(ALT == A1 & REF == A2, ALT_FREQS, (1-ALT_FREQS))) %>% # Reverse frequency if order changed!
      dplyr::mutate(
        MAF = ifelse(freq < 0.5, freq, (1-freq)), ### do we actually need MAF?
        SNP = paste0("chr", CHR, ":", BP, ":", A1, ":", A2)
      ) %>%
      dplyr::select(SNP, freq, MAF)
    df_final <- merge(df_sorted, freqs, by = "SNP")
    
    system(paste0("rm ",random.number,"*"))
  }
  
  if(!("N" %in% colnames(df_final))){
    N_hat<-median(1/((2*df_final$MAF*(1-df_final$MAF))*df_final$SE^2),na.rm = T) 
    df_final$N=ceiling(N_hat)
  }
  
# Add freq and MAF info to dataset    
  df_final <- df_final %>%
    dplyr::select("snp_original","SNP","CHR","BP","A1","A2","freq","b","varbeta","SE","P","MAF","N","type", any_of(c("s", "sdY")), "phenotype_id") %>%
    dplyr::rename(se=SE, p=P)
    
  # Remove duplicated SNPs (otherwise coloc will complain)
  setDT(df_final) # Ensure df_final is a data.table and key it by SNP
  setkey(df_final, SNP, phenotype_id)
  
  # Remove all rows with duplicated values in the SNP column, by phenotype
  df_clean <- df_final[, .SD[!duplicated(SNP) & !duplicated(SNP, fromLast =TRUE)], by = phenotype_id]
  
  return(df_clean)
}



### locus.breaker
locus.breaker <- function(
    res,
    p.sig = 5e-08,
    p.limit = 1e-05,
    hole.size = 250000,
    p.label = "P",
    chr.label = "CHR",
    pos.label = "BP"){
  
  res <- as.data.frame(res)
  res = res[order(as.numeric(res[, chr.label]), as.numeric(res[,pos.label])), ]
  res = res[which(res[, p.label] < p.limit), ]
  trait.res = c()
  
  for(j in unique(res[,chr.label])) {
    res.chr = res[which(res[, chr.label] == j), ]
    if (nrow(res.chr) > 1) {
      holes = res.chr[, pos.label][-1] - res.chr[, pos.label][-length(res.chr[,pos.label])]
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
          if (min(res.loc[, p.label]) < p.sig) {
            start.pos = min(res.loc[, pos.label], na.rm = T)
            end.pos = max(res.loc[, pos.label], na.rm = T)
            chr = j
            best.snp = res.loc[which.min(res.loc[, p.label]), 
            ]
            line.res = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (min(res.loc[, p.label]) < p.sig) {
          start.pos = min(res.loc[, pos.label], na.rm = T)
          end.pos = max(res.loc[, pos.label], na.rm = T)
          chr = j
          best.snp = res.loc[which.min(res.loc[, p.label]), 
          ]
          line.res = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (min(res.loc[, p.label]) < p.sig) {
        start.pos = min(res.loc[, pos.label], na.rm = T)
        end.pos = max(res.loc[, pos.label], na.rm = T)
        chr = j
        best.snp = res.loc[which.min(res.loc[, p.label]), 
        ]
        line.res = c(chr, start.pos, end.pos, unlist(best.snp))
        trait.res = rbind(trait.res, line.res)
      }
    }
  }
  if(!is.null(trait.res)){
    trait.res = as.data.frame(trait.res, stringsAsFactors = FALSE)
    trait.res = trait.res[, -(which(names(trait.res) == chr.label))]
    names(trait.res)[1:3] = c("chr", "start", "end")
    rownames(trait.res) <- NULL
  }
  return(trait.res)
}


#### run_dentist #### 
# Preparation of files necessary to perform DENTIST - same files will be used for COJO!!
run_dentist <- function(D=dataset_aligned
                        , locus_chr=opt$chr
                        , locus_start=opt$start
                        , locus_end=opt$end
                        , bfile="ld_panel"
                        , maf.thresh=1e-4
                        , random.number="ZUlGe4EnYqGkubYrApHu"
                        , dentist.bin="DENTIST"
){
  
  
  # Save list of snps included in the locus    
  locus_only.snp <- D %>% 
    dplyr::filter(CHR==locus_chr, BP >= locus_start, BP <= locus_end) %>%
    dplyr::pull(SNP)
  write(locus_only.snp, ncol=1,file=paste0(random.number,"_locus_only.snp.list"))
  
  # Prepare subset of plink LD files    
  exit_status = system(paste0("plink2 --bfile ", bfile," --extract ",random.number,"_locus_only.snp.list --maf ", maf.thresh, " --make-bed --out ", random.number))
  
  # Raise an error if the external command fails
  if (exit_status != 0) {
    cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
    quit(status = 1, save = "no")
  }

  #### Check if any SNP was left after extracting and filtering!
  snsp_extracted <- system(paste0("grep -E '[0-9]+ variants? remaining after main filters\\.' ", random.number, ".log"), intern = TRUE)
  snsp_extracted <- as.numeric(gsub("(\\d+) variants? remaining after main filters.", "\\1", snsp_extracted))
  
  
  if(length(snsp_extracted) > 0){
    
    # Format gwas sum stat input
    D <- D %>%
      dplyr::select("SNP","A1","A2","freq","b","se","p","N","snp_original","type", any_of(c("sdY","s")))
    
    fwrite(D %>% dplyr::select(-snp_original, -type, -any_of(c("sdY", "s"))), # to match with input required by Dentist
           file=paste0(random.number,"_sum.txt"), row.names=F,quote=F,sep="\t", na=NA)
    
    # We don't catch the exit status of the system call here, as we want to continue even if DENTIST fails
    system(paste0(dentist.bin, "/DENTIST_1.3.0.0 --gwas-summary ", random.number,"_sum.txt --bfile ", random.number, " --chrID ", locus_chr,  " --extract ", random.number, "_locus_only.snp.list --out ", random.number, " --thread-num 1"))
    
    if (file.exists(paste0(random.number, ".DENTIST.short.txt"))){ ### check that output was produced
      # Remove SNPs pointed out by DENTIST and proceed with COJO
      dentist_exclude <- fread(paste0(random.number, ".DENTIST.short.txt"), data.table = F, header = F) 
      if (nrow(dentist_exclude)>0){ ### check that output produced isn't empty
        locus_only.snp <- setdiff(locus_only.snp, dentist_exclude[,1])
        # SAVE UPDATED LIST OF SNPS!!!        
        write(locus_only.snp, ncol=1,file=paste0(random.number,"_locus_only.snp.list"))
      }
    }
  } else {
    cat(paste0(snsp_extracted, " variants remaining in the LD reference panel after SNPs extraction and MAF filter"))
    system(paste0("rm ", random.number, "*"))
    quit(save = "no", status = 0, runLast = FALSE)  # Exit the script gracefully
  }
}


### cojo.ht ###
### Performs --cojo-slct first to identify all independent SNPs and --cojo-cond then to condition upon identified SNPs
cojo.ht=function(D=dataset_aligned
                 , locus_chr=opt$chr
                 , locus_start=opt$start
                 , locus_end=opt$end
                 , p.thresh=1e-4
                 , bfile="ld_panel"
                 , maf.thresh=1e-4
                 , random.number="ZUlGe4EnYqGkubYrApHu"
                 , skip_dentist=opt$skip_dentist
){
  
  if (skip_dentist == TRUE){
    
    # Save list of snps included in the locus    
    locus_only.snp <- D %>% 
      dplyr::filter(CHR==locus_chr, BP >= locus_start, BP <= locus_end) %>%
      dplyr::pull(SNP)
    write(locus_only.snp, ncol=1,file=paste0(random.number,"_locus_only.snp.list"))
    
    # Prepare subset of plink LD files    
    exist_status = system(paste0("plink2 --bfile ", bfile," --extract ",random.number,"_locus_only.snp.list --maf ", maf.thresh, " --make-bed --out ", random.number))

    # Raise an error if the external command fails
    if (exit_status != 0) {
      cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
      quit(status = 1, save = "no")
    }

    # Format gwas sum stat input
    D <- D %>%
      dplyr::select("SNP","A1","A2","freq","b","se","p","N","snp_original","type", any_of(c("sdY","s")))
    
    fwrite(D, file=paste0(random.number,"_sum.txt"), row.names=F,quote=F,sep="\t", na=NA)
    
  }
  
  # step1 determine independent snps
  exit_status = system(paste0("gcta64 --bfile ", random.number, " --cojo-p ", p.thresh, " --maf ", maf.thresh, " --extract ", random.number, "_locus_only.snp.list --cojo-file ", random.number, "_sum.txt --cojo-slct --out ", random.number, "_step1"))
  
  # Raise an error if the external command fails
  if (exit_status != 0) {
    cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
    quit(status = 1, save = "no")
  }

  if(file.exists(paste0(random.number,"_step1.jma.cojo"))){
    dataset.list=list()
    ind.snp <- fread(paste0(random.number,"_step1.jma.cojo")) %>%
      left_join(
        D %>%
          dplyr::select(SNP,snp_original,A2,type,any_of(c("sdY", "s"))) %>%
          dplyr::rename(othA=A2),
        by="SNP")
    
    dataset.list$ind.snps <- data.frame(matrix(ncol = ncol(ind.snp), nrow = 0))
    colnames(dataset.list$ind.snps) <- colnames(ind.snp)
    dataset.list$results=list()
    
    if(nrow(ind.snp)>1){
      for(i in 1:nrow(ind.snp)){
        
        write(ind.snp$SNP[-i],ncol=1, file=paste0(random.number,"_independent.snp"))
        print(ind.snp$SNP[-i])
        
        exit_status = system(paste0("gcta64 --bfile ",random.number, " --maf ", maf.thresh, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
        
        # Raise an error if the external command fails
        if (exit_status != 0) {
          cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
          quit(status = 1, save = "no")
        }

        #### STOP ANALYSIS FOR THAT TOP SNP IN CASE OF COLLINEARITY
        if(!file.exists(paste0(random.number,"_step2.cma.cojo"))){
          cat(paste0("\n****WARNING: COJO has encountered a collinearty problem. Affected SNP will be removed from following analysis****\n\n"))
        } else {
          # Re-add type and sdY/s info, and map SNPs! AND ALSO NON EFFECT ALLELE
          step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
            left_join(
              D %>% 
                dplyr::select(SNP,snp_original,A2,type,any_of(c("sdY", "s"))) %>%
                dplyr::rename(othA=A2),
              by="SNP") %>%
            dplyr::mutate(cojo_snp=ind.snp$SNP[i])
          
          # Add SNPs to the ind.snps dataframe         
          dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp[i,])
          # Add conditioned gwas to the results list          
          dataset.list$results[[i]]=step2.res
          names(dataset.list$results)[i]=ind.snp$SNP[i]
          system(paste0("rm ",random.number,"_step2.cma.cojo"))
        }
      }
    } else {
      
      ### NB: COJO here is performed ONLY for formatting sakes - No need to condition if only one signal is found!!
      
      write(ind.snp$SNP,ncol=1,file=paste0(random.number,"_independent.snp"))
      exit_status = system(paste0("gcta64  --bfile ",random.number," --cojo-p ",p.thresh, " --maf ", maf.thresh, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
      
      # Raise an error if the external command fails
      if (exit_status != 0) {
        cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
        quit(status = 1, save = "no")
      }

      step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE)
      
      # If the locus is made up of only one SNP, step2.res will be an empty dataframe    
      if(nrow(step2.res)>0){
        step2.res <- step2.res %>%
          left_join(
            D %>% 
              dplyr::select(SNP,snp_original,A1,A2,type,any_of(c("sdY", "s"))) %>%
              dplyr::rename(othA=A2),
            by=c("SNP", "refA"="A1")
          )
      }
      
      
      #### Add back top SNP, removed from the data frame with the conditioning step
      step2.res <- rbind.fill(
        step2.res,
        ind.snp %>% dplyr::select(-bJ,-bJ_se,-pJ,-LD_r)
      )
      step2.res$cojo_snp <- ind.snp$SNP
      step2.res$bC <- step2.res$b
      step2.res$bC_se <- step2.res$se
      step2.res$pC <- step2.res$p
      
      dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp)
      dataset.list$results[[1]]=step2.res
      names(dataset.list$results)[1]=ind.snp$SNP[1]
    }
    # Remove results df possibly empty (in case of collinearity issue)
    dataset.list$results <- dataset.list$results %>% discard(is.null)
  }
  system(paste0("rm *",random.number,"*"))
  if(exists("dataset.list")){return(dataset.list)}
}


#### prep_susie_ld ####
# Prepare LD matrix for SUSIE
prep_susie_ld <- function(
    D=dataset_aligned,
    locus_chr=opt$chr,
    locus_start=opt$start,
    locus_end=opt$end,
    bfile=opt$bfile,
    maf.thresh=opt$maf,
    random.number="ZUlGe4EnYqGkubYrApHu",
    skip_dentist=opt$skip_dentist
){
  
  if (skip_dentist == TRUE){
    # Save list of snps included in the locus    
    locus_only.snp <- D %>% 
      dplyr::filter(CHR==locus_chr, BP >= locus_start, BP <= locus_end) %>%
      dplyr::pull(SNP)
    write(locus_only.snp, ncol=1,file=paste0(random.number,"_locus_only.snp.list"))
  }
  
  ### --export A include-alt --> creates a new fileset, after sample/variant filters have been applied - A: sample-major additive (0/1/2) coding, suitable for loading from R 
  exit_status = system(paste0("plink2 --bfile ", bfile, " --extract ", random.number, "_locus_only.snp.list --maf ", maf.thresh, " --export A include-alt --out ", random.number))
  
  # Raise an error if the external command fails
  if (exit_status != 0) {
    cat(paste0("Error: External command failed with exit code: ", exit_status, "\n"))
    quit(status = 1, save = "no")
  }

  geno <- fread(paste0(random.number, ".raw"))[,-c(1:6)] ### First 6 columns are FID, IID, PAT, MAT, SEX and PHENOTYPE
  
  # Check which SNPs have the same genotype for all samples and remove them
  not_same_geno <- which(sapply(geno, function(x) length(unique(x)) > 1))
  geno <- geno[, ..not_same_geno]
  
  # split the SNP names into rsID, effective and other alleles
  snp_info <- strsplit(colnames(geno), "_|\\(/|\\)") %>%
    Reduce(rbind,.) %>%
    data.frame

# If there's only one SNP in the plink .raw file, it will need further formatting
  if(ncol(snp_info)==1){
    snp_info <- data.frame(t(snp_info))
  }
  
  colnames(snp_info) <- c('SNP','ea','oa')
  rownames(snp_info) <- NULL
  colnames(geno) <- snp_info$SNP
  
  # check for which columns genotypes should be reverted
  index_to_revert <- which(!(snp_info$ea <= snp_info$oa))
  
  # Switch geno --> from 0 to 2 and vice-versa
  # Function to switch 0 to 2 and 2 to 0
  switch_0_2 <- function(x) {
    switched <- (x*-1)+2
    return(switched)
  }
  
  # Apply the transformation only to specified columns - if index_to_revert is empty, skip this step
  if(length(index_to_revert)>0){
    geno[, (index_to_revert) := lapply(.SD, switch_0_2), .SDcols = index_to_revert]
  }
  
  # Impute missing genotypes with mean value  
  geno <- apply(geno, 2, function(x) {x[is.na(x)] <- mean(x,na.rm=TRUE); return(x)})
  # Correlation matrix
  #ld <- cor(geno) #### NB: don't square it!!!!
  X_scaled <- scale(geno)  # Standardize columns
  ld <- crossprod(X_scaled) / (nrow(geno) - 1) # Same as cor(), but faster
  
  system(paste0("rm ", random.number, "*"))
  return(ld)
}



#### finemap.cojo

finemap.cojo <- function(D, cs_threshold=0.99){
  cojo_snp <- unique(D$cojo_snp)
# Format input  
    D <- D %>%
      dplyr::mutate(varbeta=bC_se^2) %>%
      dplyr::select("SNP","Chr","bp", "b", "bC","varbeta","n","pC","freq","type",any_of(c("sdY","s"))) %>%
      rename("snp"="SNP","chr"="Chr","position"="bp","beta"="bC","N"="n","pvalues"="pC","MAF"="freq")
  
  D_list <- as.list(na.omit(D)) ### move to list and keep unique value of "type" otherwise ANNOYING ERROR!
  D_list$type <- unique(D_list$type)
  if(D_list$type=="cc"){D_list$s <- unique(D_list$s)}else{D_list$sdY <- unique(D_list$sdY)}
  
# Finemap
  fine.res <- finemap.abf_NO_PRIOR(D_list) %>%
    left_join(D %>% dplyr::select(snp, b, beta, pvalues), by="snp") %>%
    dplyr::mutate(cojo_snp=cojo_snp) %>%
    dplyr::rename(bC=beta, bC_se=se, "lABF"="lABF.") %>%
    arrange(desc(SNP.PP)) %>% 
    mutate(cred.set = cumsum(SNP.PP)) %>%
# Add cojo_hit info, to merge with loci table later
    dplyr::select(snp, position, bC, bC_se, lABF, SNP.PP, cred.set, cojo_snp)

# Identify SNPs part of the credible set (as specified by cs_threshold)
  w <- which(fine.res$cred.set > cs_threshold)[1]
  cs <- fine.res %>% 
    mutate(is_cs=c(rep(TRUE, w), rep(FALSE, (nrow(fine.res)-w)))) %>%
    select(-cred.set)
  return(cs)
}



### plot.cojo.ht ###
plot.cojo.ht=function(cojo.ht.obj){
  
  if(nrow(cojo.ht.obj$ind.snps)>1){
    
    whole.dataset=c()
    for(i in 1:nrow(cojo.ht.obj$ind.snps)){
      
      tmp=cojo.ht.obj$results[[i]]
      tmp$signal=cojo.ht.obj$ind.snps$SNP[i]
      whole.dataset=rbind(whole.dataset,tmp)
    }
    
    p1 <- ggplot(cojo.ht.obj$results[[i]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_minimal()+
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=SNP),size=6,shape=23) +
      guides(fill=guide_legend(title="SNP"))
    
    p2 <- ggplot(whole.dataset,aes(x=bp,y=-log10(pC),color=signal)) +
      facet_grid(signal~.) +
      geom_point(alpha=0.8,size=3) +
      theme_minimal() +
      ggtitle("Conditioned results")
    
    p3 <- p1/p2 + plot_layout(heights = c(1, nrow(cojo.ht.obj$ind.snps)+0.2))
    
  } else {
    
    p3 <- ggplot(cojo.ht.obj$results[[1]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_minimal()+
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=SNP),size=6,shape=23)
  }
  (p3)
}



### hcolo.cojo.ht ###
hcolo.cojo.ht=function(df1 = conditional.dataset1,
                       df2 = conditional.dataset2,
                       p1=1e-4,
                       p2=1e-4,
                       p12=1e-5
                       ){
  
  df1 <- df1 %>% dplyr::rename("lABF.df1"="lABF")
  df2 <- df2 %>% dplyr::rename("lABF.df2"="lABF")
  
  p1 <- coloc:::adjust_prior(p1, nrow(df1), "1")
  p2 <- coloc:::adjust_prior(p2, nrow(df2), "2")
      
  merged.df <- merge(df1, df2, by = "snp")
  p12 <- coloc:::adjust_prior(p12, nrow(merged.df), "12")
    
  if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")
      
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- coloc:::logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)
      
  pp.abf <- coloc:::combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
      
  colo.res <- list(summary=results, results=merged.df, priors=c(p1=p1,p2=p2,p12=p12))
  class(colo.res) <- c("coloc_abf", class(colo.res))
      
## Save coloc summary        
  colo.sum <- data.frame(t(colo.res$summary))
      
## Save coloc result by SNP
  colo.full_res <- colo.res$results %>% dplyr::select(snp,lABF.df1,lABF.df2,SNP.PP.H4)

## Organise all in a list ( composed of summary + results)
  coloc.final <- list(summary=colo.sum, results=colo.full_res)
  return(coloc.final)
}






### locus.lister
locus.lister <- function(all_loci=NULL){

  pan_loci <- reduce(GRanges(seqnames = all_loci$Chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))) 
  pan_loci_non_reduced <- GRanges(seqnames = all_loci$Chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))
  
  overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) # find overlaps between all the loci and use it as an index of the unique non overlapping loci
  all_loci$pan_locus <- rep(0, nrow(all_loci)) # allocate the column 
  all_loci[overlapping@from,]$pan_locus <- overlapping@to  # assigning the number as index of which macro loci is overlapping 
  all_loci$pan_locus_name <- rep(0, nrow(all_loci))
  
  # assign a name refereed to the position of each pan_locus
  for(k in 1:length(unique(all_loci$pan_locus))){
    all_loci[which(all_loci$pan_locus==k), ]$pan_locus_name <- paste0(all_loci[which(all_loci$pan_locus==k), ]$Chr,'_',  min(all_loci[which(all_loci$pan_locus==k), ]$start), '_',  max(all_loci[which(all_loci$pan_locus==k), ]$end))
  }
  rownames(all_loci) <- NULL
  return(all_loci)
}



### coloc.plot ###
coloc.plot <- function(final.colocs.H4){
  
per.plot.data=c()
    
    for(i in 1:nrow(final.colocs.H4)){
      tmp=conditional.datasets[[final.colocs.H4$t1[i]]]$results[[final.colocs.H4$hit1[i]]]
      tmp$label=paste(final.colocs.H4$t1[i],final.colocs.H4$hit1[i],sep="-")
      tmp$group=final.colocs.H4$g1[i]
      
      if(any(!is.na(tmp$pC))){
        
        tmp=tmp[,c("bp","pC","label","group")]
        names(tmp)=c("bp","p","label","group")
        
      } else {
        
        tmp=tmp[,c("bp","p","label","group")]
      }  
      
      per.plot.data=rbind(per.plot.data,tmp)
      tmp=conditional.datasets[[x$t2[i]]]$results[[x$hit2[i]]]
      tmp$label=paste(x$t2[i],x$hit2[i],sep="-")
      tmp$group=x$g1[i]
      
      if(any(!is.na(tmp$pC))){
        
        tmp=tmp[,c("bp","pC","label","group")]
        names(tmp)=c("bp","p","label","group")
      } else {
        tmp=tmp[,c("bp","p","label","group")]
      }
      
      per.plot.data=rbind(per.plot.data,tmp)
      per.plot.data=per.plot.data[order(per.plot.data$group),]
      per.plot.data$label=factor(per.plot.data$label,levels=unique(per.plot.data$label))
      
      x$locus=locus
    }
    
    ### Plot
    pdf(paste0(outpath, "locus_", locus, "_colocalization_plot.pdf"),
        width=14, height=4*length(unique(per.plot.data$label)))  
    p1 <- ggplot(per.plot.data, aes(x=bp,y=-log10(p), fill=as.character(group))) +
      geom_point(shape=21, alpha=0.9, size=4) +
      geom_hline(yintercept = 0, linewidth=0.5) +
      facet_grid(label~.,scales="free") +
      theme_minimal() +
      #      scale_fill_manual(values=viridis(length(unique(per.plot.data$group)))) +
      #      scale_color_manual(values=viridis(length(unique(per.plot.data$group)))) +
      theme(strip.text.y.right = element_text(angle = 270, size=20),
            axis.line.y = element_line(),
            legend.position = "none") +
      ggtitle(paste("Locus ",locus," conditional regional plot"))
    print(p1)
    dev.off()
}


#### finemap.abf from coloc but modified to not include a prior/null in finemapping
finemap.abf_NO_PRIOR <- function(dataset) {
  
# Check all required input is provided  
  coloc::check_dataset(dataset,"")
  
# Compute lABF for each SNP  
  df <- coloc::process.dataset(d=dataset, suffix="")
  
# Scale  
  my.denom.log.abf <- coloc:::logsum(df$lABF)
  df$SNP.PP <- exp(df$lABF- my.denom.log.abf)
  
  return(df)
}


# Sodbo Sharapov, Human Technopole (c) 2024

#' Filters SuSiE output
#'
#' @description
#' This function filters the output from SuSiE and outputs Susie class object.
#' It return the object of class "susie" with filtered:
#' susie_output$sets$cs
#' susie_output$sets$purity
#' susie_output$sets$cs_index
#' susie_output$sets$coverage
#'
#' @param susie_output A SuSiE_R output (e.g. susieR::fitted_rss()).
#' @param cs_lbf_thr A numeric value representing the credible set log_e_BF threshold for filtering credible sets (default is 2).
#' @param pval A numeric vector representing P-values of the SNPs in the locus, tested in the susie_output.
#' @param cs_lbf_thr A numeric value representing the credible set logBF threshold for filtering credible sets (default is 2).
#' @param signal_pval_threshold A numeric value representing top SNP P-value threshold for filtering credible sets (default is 1e-6).
#' @param purity_mean_r2_threshold A numeric value representing the threshold for purity mean r2 QC metrics for filtering credible sets (default is 0.5).
#' @param purity_min_r2_threshold A numeric value representing the threshold for purity min r2 QC metrics for filtering credible sets (default is 0.5).
#' @param verbose A boolean value. If TRUE, the log information will be printed.
#'
#' @return A Susie object containing filtered fine-mapped credible sets.
#'
#' @details
#' This function takes the output from SuSiE, filters credible sets by purity metrics and CS lbf,
#' and returns a Susie object.
#'
#' @examples
#' \dontrun{
#' susie_filtered <- susie.cs.ht(
#'     susie_output = susie_output,
#'     pval = pval,
#'     cs_lbf_thr = 2,
#'     signal_pval_threshold = 1,
#'     purity_mean_r2_threshold = 0,
#'     purity_min_r2_threshold = 0,
#'     verbose = TRUE
#' )
#' }
#' @export
susie.cs.ht <- function(
    susie_output = NULL,
    pval = NULL,
    cs_lbf_thr = 2,
    signal_pval_threshold = 1,
    purity_mean_r2_threshold = 0.5,
    purity_min_r2_threshold = 0.5,
    verbose = FALSE
) {
  
  # Head of summary.susie()$cs:
  # cs cs_log10bf cs_avg_r2 cs_min_r2 variable
  # 2  22.803334 1.0000000 1.0000000 371
  # 3  13.595152 1.0000000 1.0000000 255
  # 8   4.890725 1.0000000 1.0000000 492
  # 2.738234 0.9994622 0.9994622 200,202
  
  cs_info <- susieR::summary.susie(susie_output)$cs
  
  # Get names of the CS that have lbf above threshold
  good_lbf_cs <- cs_info$cs[cs_info$cs_log10bf >= log10(exp(cs_lbf_thr))]
  
  # Get names of the CS that have mean purity above threshold
  good_purity_mean_cs <-  cs_info$cs[cs_info$cs_avg_r2 >= purity_mean_r2_threshold]
  
  # Get names of the CS that have min purity above threshold
  good_purity_min_cs <-  cs_info$cs[cs_info$cs_min_r2 >= purity_min_r2_threshold]
  
  if(!is.null(pval)) {
    # Get vector of indexes of SNPs per each CS
    min_pval_per_cs <- lapply(
      strsplit(cs_info$variable,split = ","),
      as.numeric)
    
    # Get vector of P-values of top SNP per each CS
    min_pval_per_cs <- sapply(min_pval_per_cs,function(x) min(pval[x]))
    
    # Get names of the CS that have index SNP with P-value below threshold
    good_pval_cs <- cs_info$cs[min_pval_per_cs <= signal_pval_threshold]
    
  } else {
    
    good_pval_cs <- cs_info$cs
    
  }
  
  # Get intersect of CS that satisfy all the thresholds
  cs_to_keep <-Reduce(intersect,list(good_lbf_cs, good_purity_mean_cs, good_purity_min_cs,good_pval_cs))
  
  susie_output_filtered <- susie_output
  
  # Print verbose information
  if(verbose){
    message(paste0(c("CS with lbf above threshold: ",good_lbf_cs),collapse = " "))
    message(paste0(c("CS with mean purity above threshold: ",good_purity_mean_cs),collapse = " "))
    message(paste0(c("CS with min purity threshold: ",good_purity_min_cs),collapse = " "))
    message(paste0(c("P-values of the index SNP for each CS: ",min_pval_per_cs),collapse = " "))
    message(paste0(c("CS with index SNP P-value below threshold: ",good_pval_cs),collapse = " "))
    message(paste0(c("CS which sutisfy all thresholds: ",cs_to_keep),collapse = " "))
  }
  
  # If there is at least one CS that passed threshold - report it, otherwise return NULL
  if(length(cs_to_keep)!=0){
    
    susie_output_filtered$sets$cs <- susie_output_filtered$sets$cs[paste0("L",cs_to_keep)]
    susie_output_filtered$sets$purity <- susie_output_filtered$sets$purity[paste0("L",cs_to_keep),]
    susie_output_filtered$sets$coverage <- susie_output_filtered$sets$coverage[which(susie_output_filtered$sets$cs_index %in% cs_to_keep)]
    susie_output_filtered$sets$cs_index <- cs_to_keep
    
    return(susie_output_filtered)
    
  } else {
    message("None of the credible sets satisfy selected thresholds")
    return(NULL)
  } 
  
}


### From lABF to PP - sent by Nicola
logbf_to_pp_ht = function(bf=all.coloc.join$ABF.tot) {
  denom=coloc:::logsum(bf)
  exp(bf  - denom)
}

get_beta_se_susie <- function(sus,L_index){
  
  se = sqrt(sus$mu2[L_index,] - (sus$mu[L_index,])^2)/sus$X_column_scale_factors
  beta = (sus$mu[L_index,]/sus$X_column_scale_factors)
  
  return(list(beta = beta, se = se))
  
}


# Helper function to run susie_rss with error handling
run_susie_w_tryCatch <- function(
    D_sub,
    susie_ld,
    L = L,
    coverage = coverage_value,
    max_iter = max_iter
) {  
  tryCatch({
    susie_rss(
      bhat = D_sub$BETA, 
      shat = D_sub$SE, 
      n = max(D_sub$N), 
      R = susie_ld, 
      var_y = unique(D_sub$S),
      estimate_residual_variance = FALSE,
      L = L,
      coverage = coverage,
      max_iter = max_iter
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    message(msg)
    if (grepl("The estimated prior variance is unreasonably large", msg)) {
      return("SKIP_TO_L1")
    }
    stop(e)  # Re-throw other errors
  })
}


# Run SuSiE with retries
run_susie_w_retries <- function(
    D_sub,
    susie_ld,
    L = L,
    coverage = coverage,
    max_iter = max_iter,
    min_coverage = min_coverage,
    skip_to_L1 = skip_to_L1
){
  
  fitted_rss <- NULL  # Initialize
  coverage_value_updated <- coverage
  
  # First run!  
  fitted_rss <- run_susie_w_tryCatch(
    D_sub,
    susie_ld,
    L = L,
    coverage = coverage,
    max_iter = max_iter
  )
  
  
  # "Estimated prior variance is unreasonably large" error - switch to COJO or jump to L=1?
  if (identical(fitted_rss, "SKIP_TO_L1")) {
#    message("Re-running fine-mapping with L=1")
    
#    fitted_rss <- run_susie_w_tryCatch(
#      D_sub,
#      susie_ld,
#      L = 1,
#      coverage = coverage,
#      max_iter = max_iter
#    )

#  # If cs are NULL, set the whole object to NULL - so it will be returned either as fully NULL or not NULL at all
#    if(is.null(fitted_rss$sets$cs)){
#      fitted_rss <- NULL
#    }
#    return(fitted_rss)  # return early - stop here function!
    write.table(msg, "failed_susie.txt", row.names = FALSE, col.names = FALSE)
#    quit(save = "no", status = 0, runLast = FALSE)  # Exit the script gracefully # would not work in a loop setting
    return(NULL)
  }


  # If result or cs is NULL (but not for L=1 cases) - retry lowering coverage, if still NULL jump to L=1
  while ( (is.null(fitted_rss) || is.null(fitted_rss$sets$cs)) && coverage_value_updated >= min_coverage) {
    coverage_value_updated <- coverage_value_updated - 0.01
    message("Re-running with reduced coverage: ", coverage_value_updated)
    
    fitted_rss <- run_susie_w_tryCatch(
      D_sub,
      susie_ld,
      L = L,
      coverage = coverage_value_updated,
      max_iter = max_iter
    )
  }

  # If still no credible sets, try L = 1 as a last resort
  if (coverage_value_updated < min_coverage) {
    message("Final attempt failed, reached minimum coverage of ", min_coverage, ". Re-running with L=1")
    fitted_rss <- run_susie_w_tryCatch(
      D_sub,
      susie_ld,
      L = 1,
      coverage = coverage,
      max_iter = max_iter
    )

    # If cs are still NULL, return NULL
    if(is.null(fitted_rss$sets$cs)){
      message("Region ", unique(D_sub$REGION), " for trait ", unique(D_sub$TRAITID), " was not fine-mapped.")
      return(fitted_rss)  # return early - stop here function!
    }
  }
    
    
  # Check if susie converged in given number of iterations - if not, either save for later of jump to L=1
  if (!(fitted_rss$converged)) {
    message(paste0("IBSS algorithm did not converge in ", max_iter, " iterations"))
      
    if(skip_to_L1){
      message("Re-running fine-mapping with L=1")
      fitted_rss <- run_susie_w_tryCatch(
        D_sub,
        susie_ld,
        L = 1,
        coverage = coverage,
        max_iter = max_iter
      )
    } else {
      message("Storing locus breaker intervals output")
        
      trait <- D_sub %>% pull(TRAITID) %>% unique()
      region <- D_sub %>% pull(REGION) %>% unique()
      
      inter_to_save <- D_sub %>%
        distinct(TRAITID, REGION, S) %>%
        mutate(
          CHR = gsub("(\\d+):(\\d+):(\\d+)", "\\1", REGION),
          START = gsub("(\\d+):(\\d+):(\\d+)", "\\2", REGION),
          END = gsub("(\\d+):(\\d+):(\\d+)", "\\3", REGION)
        ) %>%
        dplyr::select(TRAITID,CHR,START,END,S)
        
      fwrite(inter_to_save, file = paste0(trait, "_", region, "_retry_w_more_niters_susie_interval.csv"), sep=",", na=NA, quote=F)
      return(NULL)
    }
  }
    
  return(fitted_rss)
}
