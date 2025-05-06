include { SUSIE_FINEMAPPING       }  from "../../modules/local/susie_finemapping"
include { COJO_AND_FINEMAPPING    }  from "../../modules/local/cojo_and_finemapping"
include { APPEND_TO_MASTER_COLOC  }  from "../../modules/local/append_to_master_coloc"
include { APPEND_TO_IND_SNPS_TAB  }  from "../../modules/local/append_to_ind_snps_tab"
include { RDS_TO_ANNDATA  }  from "../../modules/local/rds_to_anndata"

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

    // Run COJO on failed SUSIE loci (only for specific errors!! Stored in the R script of susie)
    COJO_AND_FINEMAPPING(SUSIE_FINEMAPPING.out.failed_susie_loci, outdir_abspath)
    
    // Append all to independent SNPs table /// What if the channel is empty?? Can you check and behave accordingly?
    append_ind_snps = COJO_AND_FINEMAPPING.out.ind_snps_table
      .mix(SUSIE_FINEMAPPING.out.ind_snps_table)
      .groupTuple()
      .map{ tuple( it[0], it[1].flatten())}

    APPEND_TO_IND_SNPS_TAB(append_ind_snps)

    // Append all to coloc_info_master_table
    append_input_coloc = COJO_AND_FINEMAPPING.out.cojo_info_coloc_table
      .mix(SUSIE_FINEMAPPING.out.susie_info_coloc_table)
      .groupTuple()
      .map{ tuple( it[0], it[1].flatten())}

    APPEND_TO_MASTER_COLOC(append_input_coloc)

    // Collect all fine-map .rds files in AnnData
    all_rds = SUSIE_FINEMAPPING.out.susie_results_rds
      .collect()
      .map { files -> files.join(',') }
    
    RDS_TO_ANNDATA(all_rds)

  emit:
//    susie_results_rds = SUSIE_FINEMAPPING.out.susie_results_rds // to replace with AnnData
    finemap_anndata = RDS_TO_ANNDATA.out.finemap_anndata
    coloc_master = APPEND_TO_MASTER_COLOC.out.coloc_master
    ind_snps_table = APPEND_TO_IND_SNPS_TAB.out.ind_snps_table
}