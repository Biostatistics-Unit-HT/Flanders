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

    // An absent optional output produces no tasks: refinement is skipped automatically when there are no partly connected modules.
    annData_and_modules_for_cp = merged_h5ad.combine(
      IDENTIFY_REG_MODULES.out.modules_partly_connected.flatten()
    )
    COLOC_PERCOLATION(annData_and_modules_for_cp)

    // toList waits for ALL refinement tasks and emits [] if none ran.
    // Preserve actual RDS paths as process inputs, rather than hiding them in a text file containing absolute paths to other work directories.
    refinement_inputs = COLOC_PERCOLATION.out.post_cp_module_rds
      .map { result_rds, run_id_file ->
        tuple(result_rds, run_id_file.text.trim())
      }
      .toList()
      .map { records ->
        def ordered = records.sort { a, b ->
          a[0].name <=> b[0].name
        }
        def result_files = ordered.collect { it[0] }
        // Cached refinement tasks retain their IDs. A task that really re-executes writes a fresh ID, invalidating all descendants.
        def refinement_run_id = ordered
          ? ordered.collect { it[1] }.join(',')
          : 'no-refinement'
        tuple(result_files, refinement_run_id)
      }

    MASTER_TABLE(
      merged_h5ad,
      IDENTIFY_REG_MODULES.out.modules_fully_connected,
      refinement_inputs
    )

    // Carry the refinement ID through every chunk, so a fresh refinement execution invalidates joint fine-mapping even if the TSV is unchanged.
    master_table_chunks = MASTER_TABLE.out.master_table_chunk
      .flatMap { refinement_run_id, chunks ->
        def chunk_files = chunks instanceof List ? chunks : [chunks]
        chunk_files.collect { tuple(it, refinement_run_id) }
      }
    JOINT_FINEMAP_BY_CHUNK(merged_h5ad.combine(master_table_chunks))

    // Collect actual chunk RDS files, and preserve the shared refinement ID.
    joint_results = JOINT_FINEMAP_BY_CHUNK.out.joint_finemapping_rds_chunk
      .collect(flat: false)
      .map { records ->
        def ordered = records.sort { a, b -> a[1].name <=> b[1].name }
        def run_ids = ordered.collect { it[0] }.unique()
        if (run_ids.size() != 1) {
          throw new IllegalStateException('Joint fine-mapping chunks have inconsistent refinement IDs')
        }
        tuple(ordered.collect { it[1] }, run_ids[0])
      }
    FINAL_RESULTS(merged_h5ad.combine(joint_results))
}
