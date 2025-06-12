include { FIND_CS_OVERLAP_BY_CHR  } from "../../modules/local/find_cs_overlap_by_chr"
include { COLOC                   } from "../../modules/local/coloc"
include { CONCAT_ANNDATA          } from "../../modules/local/concat_anndata"
include { MAKE_COLOC_GUIDE_TABLE  } from "../../modules/local/make_coloc_guide_table"

// include { IDENTIFY_REG_MODULES     } from "../../modules/local/identify_reg_modules"

workflow RUN_COLOCALIZATION {
  take:
    credible_sets_h5ads // input channel for credible sets
		studies_to_exclude // a table containing study_id, phenotype_id to exclude from colocalization analysis

  main:
    // Concat all h5ad
    CONCAT_ANNDATA(credible_sets_h5ads)
    merged_h5ad = CONCAT_ANNDATA.out.full_anndata

    // Make a guide table, eventually filtering out previous studies
    MAKE_COLOC_GUIDE_TABLE(merged_h5ad, previous_h5ad_studies, params.coloc_filter_previous_studies)
    coloc_guide_table = MAKE_COLOC_GUIDE_TABLE.out.coloc_guide_table

    // Run COLOC process on coloc_pairs_by_batches channel
    COLOC(coloc_pairs_by_batches)

    // Collect all tables (coloc performed in batches of n pairwise tests)
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
        seed: "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tt1_study_id\tt1_phenotype_id\tt1\thit1\tt2_study_id\tt2_phenotype_id\tt2\thit2"
        ) { it.values().toList().join('\t') }
    
    coloc_results_subset.pph4
      .collectFile(
        name: "${params.coloc_id}_colocalization.table.H4.tsv",
        storeDir: "${params.outdir}/results/coloc",
        newLine: true,
        seed: "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tt1_study_id\tt1_phenotype_id\tt1\thit1\tt2_study_id\tt2_phenotype_id\tt2\thit2"
        ) { it.values().toList().join('\t') }

}