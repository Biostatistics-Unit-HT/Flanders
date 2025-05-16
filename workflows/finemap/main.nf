include { SUSIE_FINEMAPPING       }  from "../../modules/local/susie_finemapping"
include { APPEND_TO_MASTER_COLOC  }  from "../../modules/local/append_to_master_coloc"

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
    SUSIE_FINEMAPPING.out.switched_to_L1_finemap_loci_report
    .collectFile(
      keepHeader: true,
      name: "switched_to_L1_loci_report.tsv",
      storeDir: "${params.outdir}/results/switched_to_L1_loci")

    SUSIE_FINEMAPPING.out.no_credible_sets_found_loci_report
    .collectFile(
      keepHeader: true,
      name: "no_credible_sets_found_loci_report.tsv",
      storeDir: "${params.outdir}/results/not_finemapped_loci")
    
    // Append all to coloc_info_master_table
    append_input_coloc = COJO_AND_FINEMAPPING.out.cojo_info_coloc_table
      .mix(SUSIE_FINEMAPPING.out.susie_info_coloc_table)
      .groupTuple()
      .map{ tuple( it[0], it[1].flatten())}

    APPEND_TO_MASTER_COLOC(append_input_coloc)

  emit:
    susie_results_rds = SUSIE_FINEMAPPING.out.susie_results_rds
    coloc_master = APPEND_TO_MASTER_COLOC.out.coloc_master
}