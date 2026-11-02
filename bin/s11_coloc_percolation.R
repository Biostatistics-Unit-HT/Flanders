#!/usr/bin/env -S Rscript --vanilla

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(zellkonverter))
suppressMessages(library(flanders))
suppressMessages(library(hyprcoloc))
### NB: To install/load hyprcoloc library, need first to:
### module load gmp/6.3.0
### module load mpfr/4.2.1 

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--annData", default=NULL, help="AnnData cointaing fine-mapping results used for running colocalisation analysis"),
  make_option("--pph4", default=NULL, help="Posterior probability threshold for definifing H4"),
  make_option("--partly_connected_module", default=NULL, help="Colocalisation analysis output table reporting pair-wise posterior probabilities results"),
  make_option("--method", default=NULL, help="Whether to use coloc percolation or susie hypercoloc to refine partly connected modules")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# setwd("/group/pirastu/prj_014_huvec_coloc/hcoloc/modules_nf/work/ce/5559aa4615b38a72d806cd64396eff")
# opt <- list(annData="DIS_47_no_MDD_beta_fixed_2025_06_09.h5ad", pph4=0.75, partly_connected_module="partly_connected_modules_1819.tsv", method="hyprcoloc_susie")

# Read in fine-mapping anndata
ad <- zellkonverter::readH5AD(file = opt$annData, verbose = TRUE, reader = "R")

# Extract list of relevant cs
traits <- fread(opt$partly_connected_module)

# Run partially connected module through coloc percolation/hyper coloc susie to refine them

if (opt$method == "coloc_percolation"){
  
  groups <- run_coloc_percolation(
    traits,
    ad,
    pp.H4.threshold = opt$pph4,
    icoloc.reg = TRUE)

} else if (opt$method == "hyprcoloc_susie"){
 
 # Set up traits info as requested by hyercoloc function
  community <- list(
    nodes = traits$cs_name,
    chr = unique(sub(".*::(chr(?:\\d+|X|Y|M)):.+$", "\\1", traits$cs_name))
  )

  # Extract only needed chromosome from ad
  ad_chr_sub <- ad[, colData(ad)$chr == community$chr]
  ad_chr_sub <- ad_chr_sub[rowData(ad_chr_sub)$chr == community$chr, ]
  sce_per_chr <- list(ad_chr_sub)
  names(sce_per_chr) <- community$chr
  
# Run susie hypercoloc
  groups <- community2hypercoloc(sce_per_chr, community)

# Fix output!!!

### List of one element: groups$results
#### Data frame of following columns:
#[1] "iteration"                  "traits"                    
#[3] "posterior_prob"             "regional_prob"             
#[5] "candidate_snp"              "posterior_explained_by_snp"
#[7] "dropped_trait"


}

# Save coloc percolation output
old_module <- unique(traits$module)
saveRDS(groups, paste0("post_coloc_percolation_ex_", old_module, "_modules.rds"))
