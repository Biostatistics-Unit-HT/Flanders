include { IDENTIFY_REG_MODULES    } from "../../modules/local/identify_reg_modules"
include { COLOC_PERCOLATION       } from "../../modules/local/coloc_percolation"
include { MASTER_TABLE            } from "../../modules/local/master_table"
include { JOINT_FINEMAP_BY_CHUNK  } from "../../modules/local/joint_finemap_by_chunk"
include { FINAL_RESULTS           } from "../../modules/local/final_results"

workflow RUN_JOINT_FINEMAPPING {
  take:
    coloc_results_all // complete pair-wise colocalization results
    merged_h5ad       // fine-mapping AnnData used for colocalization

  main:
    IDENTIFY_REG_MODULES(coloc_results_all)

    // Emit one input pair per partly connected module. If the optional output is absent, this channel is empty and refinement is skipped.
    annData_and_modules_for_cp = merged_h5ad.combine(
      IDENTIFY_REG_MODULES.out.modules_partly_connected.flatten()
    )
    COLOC_PERCOLATION(annData_and_modules_for_cp)

    // Wait for ALL refinement tasks. toList also emits [] when none ran, allowing MASTER_TABLE to proceed without refinement results.
    // Keep actual RDS paths so Nextflow can stage and track the files.
    refinement_results = COLOC_PERCOLATION.out.post_cp_module_rds
      .flatten()
      .toList()
//      .map { files -> files.sort { a, b -> a.name <=> b.name } } // sorts the collected files alphabetically by filename

    MASTER_TABLE(
      merged_h5ad,
      IDENTIFY_REG_MODULES.out.modules_fully_connected,
      refinement_results
    )

    // Emit one task per master-table chunk.
    annData_and_master_table_by_chunk = merged_h5ad.combine(
      MASTER_TABLE.out.master_table_chunk.flatten()
    )
    JOINT_FINEMAP_BY_CHUNK(annData_and_master_table_by_chunk)

    // Wait for all joint fine-mapping results and pass the actual RDS files.
    joint_results = JOINT_FINEMAP_BY_CHUNK.out.joint_finemapping_rds_chunk
      .flatten()
      .collect()
//      .map { files -> files.sort { a, b -> a.name <=> b.name } }

    FINAL_RESULTS(merged_h5ad.combine(joint_results))
}
