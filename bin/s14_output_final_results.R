#!/usr/bin/env -S Rscript --vanilla

# Build one module-level AnnData from a pre-joint AnnData and joint-fit chunks.
# Rows remain SNPs; columns become modules. Single-CS modules keep their original
# X values. Multi-CS modules use the joint log-ABFs. Only X is exported: the
# original per-CS beta/se assays do not describe a joint module-level estimate.

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(data.table)
  library(Matrix)
})
readH5AD <- zellkonverter::readH5AD
writeH5AD <- zellkonverter::writeH5AD

# ---- Small helpers: identifiers, validation and metadata summaries ---------

module_id <- function(x) paste0("M_", sub("^M_", "", as.character(x)))

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(label, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

# Accept either historical coefficient naming convention, but reject ambiguity.
standardise_coefficients <- function(x) {
  for (old in c("bC", "bC_se")) {
    new <- c(bC = "beta", bC_se = "se")[[old]]
    if (old %in% names(x)) {
      if (new %in% names(x)) stop("Both ", old, " and ", new, " are present.")
      setnames(x, old, new)
    }
  }
  x
}

collapse_values <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_character_)
  paste(as.character(x), collapse = ", ")
}

one_value <- function(x) {
  values <- unique(as.character(x[!is.na(x)]))
  if (length(values) == 1L) values else NA_character_
}

safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

# Stable log(sum(exp(x))) without relying on coloc:::logsum.
log_sum_exp <- function(x) {
  if (anyNA(x)) stop("Missing values in a joint log-ABF row.")
  if (!length(x)) return(-Inf)
  largest <- max(x)
  if (is.infinite(largest)) return(largest)
  largest + log(sum(exp(x - largest)))
}


# Get arguments specified in the sbatch
option_list <- list(
  make_option("--all_joint_finemapping_obj", default=NULL, help="Text file containing paths and filenames of all .rds joint fine-mapping chunk results"),
  make_option("--annData", default=NULL, help="AnnData cointaing fine-mapping results used for running colocalisation analysis"),
  make_option("--coloc_id", default=NULL, help="Id code to univocally identify the colocalisation analysis"),
  make_option("--output", default=NULL, help="Output H5AD filename; default: <coloc_id>_modules_JF.h5ad")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$output)) opt$output <- paste0(opt$coloc_id, "_modules_JF.h5ad")


# Read in joint fine-mapping chunk results
all_joint_finemapping_obj_list <- readLines(opt$all_joint_finemapping_obj)
objs <- lapply(all_joint_finemapping_obj_list, readRDS) # Load all objects
message("All .rds objects loaded into R")
invisible(gc())

ad <- readH5AD(opt$annData, reader = "R", verbose = TRUE)
snps <- as.character(rowData(ad)$snp)


# ---- Define module membership once ---------------------------------------

mt <- standardise_coefficients(rbindlist(lapply(objs, `[[`, "master_table_chunk"), fill = TRUE))
require_columns(mt, c("module", "cs_name", "chr", "start", "end", "phenotype_id", "study_id"), "Master table")
mt[, module := module_id(module)]

sizes <- mt[, .(n_cs = .N), by = module]
single_ids <- sizes[n_cs == 1L, module]
multi_ids <- sizes[n_cs > 1L, module]
final_ids <- c(single_ids, multi_ids)

# Optional metadata have explicit missing values instead of get0()/scope lookup.
for (field in c("most_likely_snp", "most_likely_beta", "most_likely_se")) {
  if (!field %in% names(mt)) mt[, (field) := NA]
}

multi_cd <- mt[module %in% multi_ids, .(
  cs_name = module[1L],
  list_of_cs = collapse_values(cs_name),
  chr = one_value(chr),
  start = safe_min(start),
  end = safe_max(end),
  phenotype_id = collapse_values(unique(phenotype_id)),
  most_likely_snp = one_value(most_likely_snp),
  list_of_most_likely_beta = collapse_values(most_likely_beta),
  list_of_most_likely_se = collapse_values(most_likely_se),
  study_id = one_value(study_id),
  list_of_study_id = collapse_values(unique(study_id)),
  n_cs = .N
), by = module]
multi_cd <- multi_cd[match(multi_ids, module)]


# ---- Single modules: copy original CS columns and retain their metadata ----

# Subset by the membership table, rather than unique(cs_name).
single_map <- mt[match(single_ids, module), .(module, cs_name)]
single_idx <- match(single_map$cs_name, colnames(ad))
if (anyNA(single_idx)) stop("Some singleton CS identifiers are absent from AnnData.")
single_x <- assay(ad, "X")[, single_idx, drop = FALSE]
colnames(single_x) <- single_ids
single_cd <- standardise_coefficients(as.data.table(
  as.data.frame(colData(ad)[single_idx, , drop = FALSE])
))
single_fields <- c("study_id", "snp", "beta", "se")
if (length(single_ids)) {
  require_columns(single_cd, single_fields, "Single-CS metadata")
} else {
  # No singleton data are needed, but retain an empty table with valid fields.
  for (field in setdiff(single_fields, names(single_cd))) {
    single_cd[, (field) := character()]
  }
}
single_cd[, `:=`(
  module = single_ids,
  cs_name = single_ids,
  list_of_cs = single_map$cs_name,
  list_of_study_id = as.character(study_id),
  n_cs = 1L,
  most_likely_snp = as.character(snp),
  list_of_most_likely_beta = as.character(beta),
  list_of_most_likely_se = as.character(se)
)]
# Preserve the original deliberate removal of CS-specific L/type annotations.
remove_fields <- intersect(c("L", "type"), names(single_cd))
if (length(remove_fields)) single_cd[, (remove_fields) := NULL]

# ---- Multi modules: map sparse chunk entries onto the shared SNP axis ------

triplets <- vector("list", length(objs))
joint_meta <- vector("list", length(objs))
seen_modules <- character()

for (k in if (length(multi_ids)) seq_along(objs) else integer()) {
  m <- objs[[k]]$joint.lABF.matrix
  if (!is(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  rownames(m) <- module_id(rownames(m))
  selected <- which(rownames(m) %in% multi_ids)
  if (!length(selected)) next
  m <- m[selected, , drop = FALSE]
  
  seen_modules <- c(seen_modules, rownames(m))
  snp_idx <- match(colnames(m), snps)
  entries <- summary(m)
  triplets[[k]] <- data.table(
    i = snp_idx[entries$j],
    j = match(rownames(m), multi_ids)[entries$i], x = entries$x
  )
  
  meta <- as.data.table(copy(objs[[k]]$cs_metadata))
  meta[, module := module_id(module)]
  idx <- match(rownames(m), meta$module)

  joint_meta[[k]] <- data.table(
    module = rownames(m), min_res_labf = meta$min.labf[idx],
    logsum.logABF = vapply(seq_len(nrow(m)), function(j) {
      log_sum_exp(as.numeric(m[j, , drop = FALSE]))
    }, numeric(1))
  )
  message("Processed chunk ", k, "/", length(objs))
}

if (length(multi_ids)) {
  entries <- rbindlist(triplets)
  multi_x <- sparseMatrix(i = entries$i, j = entries$j, x = entries$x,
                          dims = c(length(snps), length(multi_ids)),
                          dimnames = list(rownames(ad), multi_ids))
  meta <- rbindlist(joint_meta)
  idx <- match(multi_ids, meta$module)
  multi_cd[, `:=`(min_res_labf = meta$min_res_labf[idx],
                  logsum.logABF = meta$logsum.logABF[idx], snp = most_likely_snp)]
} else {
  multi_x <- sparseMatrix(i = integer(), j = integer(), x = numeric(),
                          dims = c(length(snps), 0L))
}

# ---- Assemble one SCE with one consistent module-metadata table ------------

# Keep singleton-only annotations, filling absent multi-module fields with NA.
module_cd <- rbindlist(list(single_cd, multi_cd), use.names = TRUE, fill = TRUE)
module_cd <- module_cd[match(final_ids, module)]
module_cd[, list_of_phenotype_id := as.character(phenotype_id)] # phenotype_id already contains the combined IDs for multi-CS modules

# SNP IDs follow chr:pos:effect_allele:other_allele. For multi-CS modules, use the representative most_likely_snp.
multi_rows <- which(module_cd$cs_n %in% multi_ids)

if (length(multi_rows)) {
  alleles <- tstrsplit(
    as.character(module_cd$most_likely_snp[multi_rows]),
    ":",
    fixed = TRUE,
    keep = c(3L, 4L)
  )
  
  module_cd[multi_rows, `:=`(
    a1 = alleles[[1L]],
    a0 = alleles[[2L]]
  )]
}


# Remove redundant fields from the final module metadata
remove_fields <- intersect(
  c("comment_section", "beta", "se", "study_id", "phenotype_id", "module"),
  names(module_cd)
)

if (length(remove_fields)) {
  module_cd[, (remove_fields) := NULL]
}


module_cd <- as.data.frame(module_cd)
rownames(module_cd) <- final_ids
final_x <- cbind(single_x, multi_x)
colnames(final_x) <- final_ids
rownames(final_x) <- rownames(ad)
sce_final <- SingleCellExperiment(
  assays = list(X = final_x), rowData = rowData(ad),
  colData = S4Vectors::DataFrame(module_cd)
)

stopifnot(identical(colnames(sce_final), rownames(colData(sce_final))),
          identical(rownames(sce_final), rownames(ad)),
          ncol(sce_final) == nrow(sizes))
validObject(sce_final)


# ---- Outputs  ----------------------------------

fwrite(mt, paste0(opt$coloc_id, "_master_table.mt_all.tsv"), sep = "\t", quote = FALSE, na = "NA")
writeH5AD(sce_final, file = opt$output, X_name = "X")

message("Saved ", opt$output, "\n", nrow(sce_final), " SNPs; ", length(single_ids), " single-CS and ", length(multi_ids), " multi-CS modules")
