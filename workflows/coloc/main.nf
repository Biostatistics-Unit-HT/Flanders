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
        name: "${params.coloc_id}_colocalization.table.all.tsv",
        storeDir: "${params.outdir}/results/coloc",
        keepHeader: true, skip: 1)
      // .combine(credible_sets)
    

    // Split the coloc_results_all channel into two branches based on the pph4 and pph3 thresholds
    // The resulting subsets are also saved to tables
    coloc_results_all
      .splitCsv(header:true, sep:"\t")
      .branch { row ->
        pph4: row['PP.H4.abf'].toFloat() >= params.pph4_threshold
        pph3: row['PP.H3.abf'].toFloat() >= params.pph3_threshold
      }
      .set { coloc_results_subset }

    coloc_results_subset.pph3
      .collectFile(
        name: "${params.coloc_id}_colocalization.table.H3.tsv",
        storeDir: "${params.outdir}/results/coloc",
        newLine: true,
        seed: "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tt1_study_id\tt1\thit1\tt2_study_id\tt2\thit2"
        ) { it.values().toList().join('\t') }
    
    coloc_results_subset.pph4
      .collectFile(
        name: "${params.coloc_id}_colocalization.table.H4.tsv",
        storeDir: "${params.outdir}/results/coloc",
        newLine: true,
        seed: "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tt1_study_id\tt1\thit1\tt2_study_id\tt2\thit2"
        ) { it.values().toList().join('\t') }

    // Run COLOC process on coloc_pairs_by_batches channel
    // IDENTIFY_REG_MODULES(coloc_results_all)
}