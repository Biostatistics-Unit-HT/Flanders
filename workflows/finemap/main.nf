include { SUSIE_FINEMAPPING       }  from "../../modules/local/susie_finemapping"
include { APPEND_TO_MASTER_COLOC  }  from "../../modules/local/append_to_master_coloc"
include { RDS_TO_ANNDATA          }  from "../../modules/local/rds_to_anndata"

workflow RUN_FINEMAPPING {
  take:
    finemap_configuration // configuration for each study finemap
    finemapped_loci // finemapped loci channel from MUNG_AND_LOCUS_BREAKER
    munged_stats // output channel of MUNG_AND_LOCUS_BREAKER
    outdir_abspath // value with absolute path to output directory

  main:
    finemapping_input = finemap_configuration
    .combine(finemapped_loci, by:0)
    .combine(munged_stats, by:0)

    // Run SUSIE_FINEMAPPING process on finemapping_input channel
    SUSIE_FINEMAPPING(finemapping_input, outdir_abspath)

    // Concatenate all susie switched to L=1 and no credible sets found fine-mapping loci and publish it
    SUSIE_FINEMAPPING.out.finemapped_L1_prior_variance_too_large
    .collectFile(
      keepHeader: true,
      name: "FINEMAPPED_L1_prior_variance_too_large.tsv",
      storeDir: "${params.outdir}/results/finemapping_exceptions")

    SUSIE_FINEMAPPING.out.finemapped_L1_IBSS_algorithm_did_not_converge
    .collectFile(
      keepHeader: true,
      name: "FINEMAPPED_L1_IBSS_algorithm_did_not_converge.tsv",
      storeDir: "${params.outdir}/results/finemapping_exceptions")

    SUSIE_FINEMAPPING.out.not_finemapped_no_credible_sets_found
    .collectFile(
      keepHeader: true,
      name: "NOT_FINEMAPPED_no_credible_sets_found.tsv",
      storeDir: "${params.outdir}/results/finemapping_exceptions")
    
    // Append all to coloc_info_master_table
    append_input_coloc = SUSIE_FINEMAPPING.out.susie_info_coloc_table
      .groupTuple()
      .map{ tuple( it[0], it[1].flatten())}

    APPEND_TO_MASTER_COLOC(append_input_coloc)
    
    // Collect all fine-map .rds files in AnnData
    all_rds = SUSIE_FINEMAPPING.out.susie_results_rds
      .collect()
      
    RDS_TO_ANNDATA(all_rds)

  emit:
    finemap_anndata = RDS_TO_ANNDATA.out.finemap_anndata
    susie_results_rds = SUSIE_FINEMAPPING.out.susie_results_rds
    coloc_master = APPEND_TO_MASTER_COLOC.out.coloc_master
}