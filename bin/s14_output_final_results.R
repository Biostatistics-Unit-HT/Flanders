#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(anndata))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--all_joint_finemapping_obj", default=NULL, help="List of paths and filenames of all .rds joint fine-mapping chunk results"),
  make_option("--annData", default=NULL, help="AnnData cointaing fine-mapping results used for running colocalisation analysis"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# -------------------------
# Load inputs
# -------------------------
# Read in joint fine-mapping chunk results
all_joint_finemapping_obj_list <- readLines(opt$all_joint_finemapping_obj)
objs <- lapply(all_joint_finemapping_obj_list, readRDS) # Load all objects
message("All .rds objects loaded into R")
invisible(gc())

for (i in seq_along(objs)) {
  mat <- objs[[i]]$joint.lABF.matrix
  logsums <- vapply(seq_len(nrow(mat)),
                    \(j) coloc:::logsum(mat[j,]),
                    numeric(1))
  objs[[i]]$cs_metadata$logsum.logABF <- logsums
  if (i %% 10 == 0) message("Processed chunk ", i, "/", length(objs))
}

# Pre-joint SCE (columns are CS)
ad <- zellkonverter::readH5AD(file = opt$annData, verbose = TRUE, reader = "R")
ad_vars <- SingleCellExperiment::rowData(ad)  # SNP metadata for rows of the module-level SCE


# -------------------------
# Build full master table; identify single vs multi modules
# -------------------------
mt_all <- data.table::rbindlist(lapply(objs, `[[`, "master_table_chunk"), fill = TRUE)
mt_all[, module := as.character(module)]

fwrite(mt_all, file = paste0(opt$coloc_id, "_master_table.mt_all.tsv"), quote=F, na=NA, sep="\t")

module_sizes <- mt_all[, .N, by = module]
single_mods  <- module_sizes[N == 1L, module]
multi_mods   <- module_sizes[N  > 1L, module]

# CS lists for those module sets
cs_single <- mt_all[module %in% single_mods, unique(cs_name)]
cs_multi  <- mt_all[module %in% multi_mods,  unique(cs_name)]

# -------------------------
# Create pre-joint SCE with ONLY single-module CS
# -------------------------
stopifnot(length(cs_single) > 0)
missing_single <- setdiff(cs_single, colnames(ad))
if (length(missing_single)) {
  stop(sprintf("These single-CS are not in colnames(ad_cs): %s",
               paste(missing_single, collapse = ", ")))
}
sce_single <- ad[, cs_single, drop = FALSE]

# -------------------------
# Build joint matrix ONLY for multi-CS modules
# -------------------------
# Module names for joint matrix are ONLY multi_mods (order fixed by unique())
module_names <- unique(mt_all[module %in% multi_mods, module])
module_lookup <- setNames(seq_along(module_names), module_names)

# SNP universe/order comes from the pre-joint AnnData (keeps consistency)
snp_universe <- as.character(ad_vars$snp)
snp_lookup   <- setNames(seq_along(snp_universe), snp_universe)

# Two-pass style: pre-size to total nnz; we will trim after filtering
nnz_per_obj <- vapply(objs, function(o) {
  M <- o$joint.lABF.matrix
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  M@p[length(M@p)]
}, integer(1))
nnz_total <- sum(nnz_per_obj)

I <- integer(nnz_total)
J <- integer(nnz_total)
V <- numeric(nnz_total)
ptr <- 1L

for (k in seq_along(objs)) {
  print(k)
  M <- objs[[k]]$joint.lABF.matrix
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  
  # local module rows and which of them are multi-mods
  rnames <- rownames(M)
  # map only rows that belong to multi_mods; others become NA and will be dropped
  local2global <- unname(as.integer(module_lookup[rnames]))
  
  # sparse slots
  p  <- M@p
  ii <- M@i
  xx <- M@x
  
  # repeat column indices by nnz in each column
  col_rep <- rep.int(seq_len(ncol(M)), diff(p))
  
  # map SNP columns to global indices (NA if SNP not in universe)
  snp_idx_by_col <- unname(snp_lookup[colnames(M)])
  Jk <- snp_idx_by_col[col_rep]
  
  # map nonzero row indices (0-based) to global module indices
  Ik <- local2global[ii + 1L]
  
  keep <- !is.na(Ik) & !is.na(Jk)
  nk <- sum(keep)
  if (nk > 0L) {
    rng <- ptr:(ptr + nk - 1L)
    I[rng] <- Ik[keep]
    J[rng] <- Jk[keep]
    V[rng] <- xx[keep]
    ptr <- ptr + nk
  }
}

# Trim to the nnz actually kept (multi modules + in-universe SNPs)
if (ptr <= nnz_total) {
  I <- I[seq_len(ptr - 1L)]
  J <- J[seq_len(ptr - 1L)]
  V <- V[seq_len(ptr - 1L)]
}

# Construct sparse matrix (modules x SNPs), ONLY multi-CS modules
lABF_matrix_sparse <- sparseMatrix(
  i = I, j = J, x = V,
  dims = c(length(module_names), length(snp_universe)),
  dimnames = list(module_names, snp_universe)
)

lABF_matrix_sparse <- t(lABF_matrix_sparse)

colnames(lABF_matrix_sparse) <- paste0("M_",colnames(lABF_matrix_sparse))

gc()

# -------------------------
# Build module-level SCE (rows = SNPs, cols = modules) for MULTI modules
# -------------------------
# Summarize one colData row per module, restricting to multi_mods
summ_coldata <- mt_all[module %in% multi_mods,
                       {
                         list_of_cs               <- paste(cs_name, collapse = ", ")
                         chr_val                  <- { u <- unique(as.character(chr)); if (length(u) == 1L) u else NA_character_ }
                         start_val                <- suppressWarnings(min(start, na.rm = TRUE))
                         end_val                  <- suppressWarnings(max(end,   na.rm = TRUE))
                         phenotype_id_val         <- paste(unique(as.character(phenotype_id)), collapse = ", ")
                         mls_u                    <- unique(as.character(get0("most_likely_snp", ifnotfound = NA)))
                         most_likely_snp_val      <- if (length(mls_u) == 1L) mls_u else NA_character_
                         list_of_most_likely_beta <- if ("most_likely_beta" %in% names(.SD)) paste(most_likely_beta, collapse = ", ") else NA_character_
                         list_of_most_likely_se   <- if ("most_likely_se"   %in% names(.SD)) paste(most_likely_se,   collapse = ", ") else NA_character_
                         study_ids                <- unique(as.character(study_id))
                         study_id_val             <- if (length(study_ids) == 1L) study_ids else NA_character_
                         list_of_study_id         <- paste(study_ids, collapse = ", ")
                         
                         data.table(
                           # DO NOT add `module = ...` here; it’s already provided by `by = module`
                           list_of_cs               = list_of_cs,
                           chr                      = chr_val,
                           start                    = start_val,
                           end                      = end_val,
                           phenotype_id             = phenotype_id_val,
                           most_likely_snp          = most_likely_snp_val,
                           list_of_most_likely_beta = list_of_most_likely_beta,
                           list_of_most_likely_se   = list_of_most_likely_se,
                           study_id                 = study_id_val,
                           list_of_study_id         = list_of_study_id,
                           n_cs                     = .N
                         )
                       },
                       by = module
]

col_df <- as.data.frame(summ_coldata)
col_df$module <- paste0("M_",col_df$module)
rownames(col_df) <- col_df$module
# align to module_names order used in lABF_matrix_sparse
col_df <- col_df[paste0("M_",module_names), , drop = FALSE]

sce_multi_modules <- SingleCellExperiment(
  assays  = list(X = lABF_matrix_sparse),  # rows = SNPs, cols = MULTI modules
  rowData = as.data.frame(ad_vars),
  colData = col_df
)

assays(sce_multi_modules)[["beta"]] <- Matrix(
  0, nrow = nrow(sce_multi_modules), ncol = ncol(sce_multi_modules),
  dimnames = dimnames(assay(sce_multi_modules, "X")),
  sparse = TRUE
)

assays(sce_multi_modules)[["se"]] <- Matrix(
  0, nrow = nrow(sce_multi_modules), ncol = ncol(sce_multi_modules),
  dimnames = dimnames(assay(sce_multi_modules, "X")),
  sparse = TRUE
)

# -------------------------
# Combine two sce
# -------------------------

## --- Fix sce_single ---
cds <- as.data.frame(colData(sce_single), stringsAsFactors = FALSE)
cds$list_of_cs               <- cds$cs_name
cds$list_of_study_id         <- cds$study_id
cds$n_cs                     <- 1L
cds$most_likely_snp          <- cds$snp
cds$list_of_most_likely_beta <- cds$bC
cds$list_of_most_likely_se   <- cds$bC_se
cds <- cds %>% rename("beta"="bC","se"="bC_se")
cds$type <- NULL
cds$panel <- NA
cds$category <- NA
colData(sce_single) <- DataFrame(cds)

## --- Fix sce_multi_modules ---
cdm <- as.data.frame(colData(sce_multi_modules), stringsAsFactors = FALSE)
cdm$category        <- NA
cdm$top_pvalue      <- NA
cdm$panel           <- NA
cdm$cs_name         <- cdm$module
meta_all <- data.table::rbindlist(lapply(objs, `[[`, "cs_metadata"))
meta_all[, module := as.character(module)]
meta_all$module <- paste0("M_", meta_all$module)
cdm$min_res_labf <- meta_all$min.labf[match(as.character(cdm$module), meta_all$module)]
cdm$logsum.logABF <- meta_all$logsum.logABF[match(as.character(cdm$module), meta_all$module)]
cdm$snp <- cdm$most_likely_snp
parts <- tstrsplit(cdm$most_likely_snp, ":", fixed=TRUE)
cdm$a1 <- parts[[3]]; cdm$a0 <- parts[[4]]
cdm$freq <- cdm$N <- cdm$beta <- cdm$se <- cdm$min.abs.corr <- cdm$mean.abs.corr <- cdm$median.abs.corr <- NA
cdm$comment_section <- NA
cdm$module <- NULL
colData(sce_multi_modules) <- DataFrame(cdm)

## --- modify sce_multi_modules ---
# add constant coverage = 0.99
colData(sce_multi_modules)$coverage <- 0.99

## --- modify sce_single ---
# remove L
colData(sce_single)$L <- NULL

# add empty comment_section
colData(sce_single)$comment_section <- NA

# remove category
colData(sce_single)$category <- NA

sce_final <- cbind(sce_single,sce_multi_modules)

# -------------------------
# Save outputs
# -------------------------
zellkonverter::writeH5AD(
  sce_final,
  file = paste0(opt$coloc_id, "_mndules_JF.h5ad")
  )

message("Done")
