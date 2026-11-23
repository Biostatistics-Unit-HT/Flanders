include { SUSIE_FINEMAPPING                         }  from "../../modules/local/susie_finemapping"
include { APPEND_TO_MASTER_COLOC                    }  from "../../modules/local/append_to_master_coloc"
include { RDS_TO_ANNDATA                            }  from "../../modules/local/rds_to_anndata"
include { PUBLISH_FILE as PUBLISH_L1_PRIOR_VARIANCE } from "../../modules/local/publish_file"
include { PUBLISH_FILE as PUBLISH_L1_NONCONVERGENCE } from "../../modules/local/publish_file"
include { PUBLISH_FILE as PUBLISH_L1_QC_RECOVERY    } from "../../modules/local/publish_file"
include { PUBLISH_FILE as PUBLISH_NO_CREDIBLE_SETS  } from "../../modules/local/publish_file"
include { PUBLISH_FILE as PUBLISH_NO_LD_VARIANTS    } from "../../modules/local/publish_file"

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
    SUSIE_FINEMAPPING(
      finemapping_input,
      outdir_abspath,
      params.susie_max_iter,
      params.susie_qc_cs_bf_thr,
      params.susie_qc_pval_thr,
      params.susie_qc_mean_r2_thr,
      params.susie_qc_min_r2_thr,
      params.publish_susie
    )

    // Concatenate all susie switched to L=1 and no credible sets found fine-mapping loci and publish it
    l1_prior_variance_report = SUSIE_FINEMAPPING.out.finemapped_L1_prior_variance_too_large
    .collectFile(
      keepHeader: true,
      name: "FINEMAPPED_L1_prior_variance_too_large.tsv")

    l1_nonconvergence_report = SUSIE_FINEMAPPING.out.finemapped_L1_IBSS_algorithm_did_not_converge
    .collectFile(
      keepHeader: true,
      name: "FINEMAPPED_L1_IBSS_algorithm_did_not_converge.tsv")

    l1_qc_recovery_report = SUSIE_FINEMAPPING.out.finemapped_L1_recover_after_susie_QC
    .collectFile(
      keepHeader: true,
      name: "FINEMAPPED_L1_recover_after_susie_QC.tsv")

    no_credible_sets_report = SUSIE_FINEMAPPING.out.not_finemapped_no_credible_sets_found
    .collectFile(
      keepHeader: true,
      name: "NOT_FINEMAPPED_no_credible_sets_found.tsv")

    no_ld_variants_report = SUSIE_FINEMAPPING.out.not_finemapped_no_variants_from_locus_in_LD_ref
    .collectFile(
      keepHeader: true,
      name: "NOT_FINEMAPPED_no_variants_from_locus_in_LD_ref.tsv")
    
    // Append all to coloc_info_master_table
    append_input_coloc = SUSIE_FINEMAPPING.out.susie_info_coloc_table
      .groupTuple()
      .map{ tuple( it[0], it[1].flatten())}

    PUBLISH_L1_PRIOR_VARIANCE(l1_prior_variance_report, "${params.outdir}/results/finemapping_exceptions")
    PUBLISH_L1_NONCONVERGENCE(l1_nonconvergence_report, "${params.outdir}/results/finemapping_exceptions")
    PUBLISH_L1_QC_RECOVERY(l1_qc_recovery_report, "${params.outdir}/results/finemapping_exceptions")
    PUBLISH_NO_CREDIBLE_SETS(no_credible_sets_report, "${params.outdir}/results/finemapping_exceptions")
    PUBLISH_NO_LD_VARIANTS(no_ld_variants_report, "${params.outdir}/results/finemapping_exceptions")

    APPEND_TO_MASTER_COLOC(append_input_coloc)
    
    // Collect all fine-map .rds files in AnnData
    all_rds = SUSIE_FINEMAPPING.out.susie_results_rds
      .collect()
      
    RDS_TO_ANNDATA(all_rds, params.finemap_id)

  emit:
    finemap_anndata = RDS_TO_ANNDATA.out.finemap_anndata
}