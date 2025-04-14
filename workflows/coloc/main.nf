include { FIND_CS_OVERLAP_BY_CHR  } from "../../modules/local/find_cs_overlap_by_chr"
include { COLOC                   } from "../../modules/local/coloc"
// include { IDENTIFY_REG_MODULES     } from "../../modules/local/identify_reg_modules"

workflow RUN_COLOCALIZATION {
  take:
    credible_sets // input channel for credible sets
  
  main:
    // Ensure the folder to store not finemapped loci exists
    file("${params.outdir}/results/coloc").mkdirs()

    // Run FIND_CS_OVERLAP process on all_cond_datasets_cs channel
    FIND_CS_OVERLAP_BY_CHR(credible_sets)

    // Define input channel for performing pair-wise colocalisation analysis (in batches)
    coloc_pairs_by_batches = FIND_CS_OVERLAP_BY_CHR.out.coloc_pairwise_guide_table
      .splitText(by:params.coloc_batch_size, keepHeader:true, file:true)

    // Run COLOC process on coloc_pairs_by_batches channel
    COLOC(coloc_pairs_by_batches)

    // Define input channel for identifying regulatory modules - collect all tables (coloc performed in batches of n pairwise tests)
    coloc_results_all = COLOC.out.colocalization_table_all_by_chunk
      .collectFile(
        name: 'colocalization_table_all_merged.txt',
        storeDir: "${params.outdir}/results/coloc",
        keepHeader: true, skip: 1)
      .combine(credible_sets)
    
    // Run COLOC process on coloc_pairs_by_batches channel
    // IDENTIFY_REG_MODULES(coloc_results_all)
}