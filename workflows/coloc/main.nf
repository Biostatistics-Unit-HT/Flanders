include { CONCAT_ANNDATA          } from "../../modules/local/concat_anndata"
include { FIND_CS_OVERLAP_BY_CHR  } from "../../modules/local/find_cs_overlap_by_chr"
include { MAKE_COLOC_GUIDE_TABLE  } from "../../modules/local/make_coloc_guide_table"
include { COLOC                   } from "../../modules/local/coloc"
include { IDENTFY_REG_MODULES     } from "../../modules/local/identify_reg_modules"
include { COLOC_PERCOLATION       } from "../../modules/local/coloc_percolation"
include { MASTER_TABLE            } from "../../modules/local/master_table"
include { JOINT_FINEMAP_BY_CHUNK  } from "../../modules/local/joint_finemap_by_chunk"
include { FINAL_RESULTS           } from "../../modules/local/final_results"

workflow RUN_COLOCALIZATION {
  take:
    credible_sets_h5ads // a list of h5ad files with credible sets
		previous_studies_from_h5ad // a table containing study_id, phenotype_id from previous h5ad to exclude
    exclude_studies_file // a table containing study_id, phenotype_id to exclude

  main:
    // Concat all h5ad if we have more than one anndata
    credible_sets_h5ads
      .branch { files ->
          single: files.size() == 1
              return files[0]  // Return the single file
          multiple: files.size() > 1
              return files     // Return the list of files
      }
    .set { h5ad_channel }

    CONCAT_ANNDATA(h5ad_channel.multiple)
    merged_h5ad = h5ad_channel.single.mix(CONCAT_ANNDATA.out.full_anndata)

    if (params.coloc_guide_table) {
      // If guide table is provided read this and ignore guide table from anndata
      coloc_guide_table = file(params.coloc_guide_table, checkIfExist: true)
    } else {
      // Make a guide table, eventually filtering out previous studies
      MAKE_COLOC_GUIDE_TABLE(merged_h5ad, previous_studies_from_h5ad, exclude_studies_file)
      coloc_guide_table = MAKE_COLOC_GUIDE_TABLE.out.coloc_guide_table
    }

    coloc_pairs_by_batches = coloc_guide_table
        .splitText(by: params.coloc_batch_size, keepHeader: true, file: true)

    coloc_input_ch = merged_h5ad.combine(coloc_pairs_by_batches)

    // Run COLOC process on coloc_pairs_by_batches channel
    COLOC(coloc_input_ch)

    // Collect all tables (coloc performed in batches of n pairwise tests)
    coloc_results_all = COLOC.out.colocalization_chunk
      .collectFile(
        name: "${params.coloc_id}_colocalization.table.all.tsv",
        storeDir: "${params.outdir}/results/coloc",
        keepHeader: true, skip: 1)

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
        seed: "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tt1_study_id\tt1\tt2_study_id\tt2\thit1\thit2"
        ) { it.values().toList().join('\t') }
    
    coloc_results_subset.pph4
      .collectFile(
        name: "${params.coloc_id}_colocalization.table.H4.tsv",
        storeDir: "${params.outdir}/results/coloc",
        newLine: true,
        seed: "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tt1_study_id\tt1\tt2_study_id\tt2\thit1\thit2"
        ) { it.values().toList().join('\t') }

    
    // Identify co-regulated modules
    IDENTFY_REG_MODULES(coloc_results_all)

   // Handle COLOC_PERCOLATION conditionally
    modules_partly_connected = IDENTFY_REG_MODULES.out.modules_partly_connected
    
    // Create a default empty file that will be used if no COLOC_PERCOLATION runs
    empty_post_cp_file = Channel.value("# No post-coloc percolation modules")
        .collectFile(name: 'empty_post_cp_module_rds.txt', newLine: true)

    // Run COLOC_PERCOLATION only if modules exist
    annData_and_modules_for_cp = merged_h5ad.combine(modules_partly_connected.flatten())
    COLOC_PERCOLATION(annData_and_modules_for_cp)

    // Collect COLOC_PERCOLATION outputs if they exist
    post_cp_from_process = COLOC_PERCOLATION.out.post_cp_module_rds
        .map { file -> file.toString() }
        .collectFile(name: 'all_post_cp_module_rds.txt', newLine: true)

    // Use the process output if available, otherwise use the empty file
    all_post_cp_module_rds = post_cp_from_process
        .mix(empty_post_cp_file)
        .first()

    // Create the master table and perform joint fine-mapping
    MASTER_TABLE(merged_h5ad, IDENTFY_REG_MODULES.out.modules_fully_connected, all_post_cp_module_rds)

    // Split chuncked master table into individual emissions and combine with annData
    annData_and_master_table_by_chunk = merged_h5ad.combine(MASTER_TABLE.out.master_table_chunk.flatten())

    // Perform joint fine-mapping by chunk
    JOINT_FINEMAP_BY_CHUNK(annData_and_master_table_by_chunk)

    // Collect ALL .rds outputs from JOINT_FINEMAP_BY_CHUNK and write to a text file
    all_joint_finemapping_obj = JOINT_FINEMAP_BY_CHUNK.out.joint_finemapping_rds_chunk
        .map { file -> file.toString() }  // Convert path to string
        .collectFile(name: 'all_joint_finemapping_rds.txt', newLine: true)

    // Collate all .rds together
    FINAL_RESULTS(merged_h5ad.combine(all_joint_finemapping_obj))

}