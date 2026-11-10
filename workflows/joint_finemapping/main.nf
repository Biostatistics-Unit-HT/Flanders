include { IDENTIFY_REG_MODULES    } from "../../modules/local/identify_reg_modules"
include { COLOC_PERCOLATION       } from "../../modules/local/coloc_percolation"
include { MASTER_TABLE            } from "../../modules/local/master_table"
include { JOINT_FINEMAP_BY_CHUNK  } from "../../modules/local/joint_finemap_by_chunk"
include { FINAL_RESULTS           } from "../../modules/local/final_results"

workflow RUN_JOINT_FINEMAPPING {
  take:
    coloc_results_all // pairwise tests coloc table results (no PP.H3/4 filtering)
    merged_h5ad       // anndata of all credible sets submitted to colocalization

  main:
    // Identify co-regulated modules
    IDENTIFY_REG_MODULES(coloc_results_all)

   // Handle COLOC_PERCOLATION conditionally
    modules_partly_connected = IDENTIFY_REG_MODULES.out.modules_partly_connected
    
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
    MASTER_TABLE(merged_h5ad, IDENTIFY_REG_MODULES.out.modules_fully_connected, all_post_cp_module_rds)

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