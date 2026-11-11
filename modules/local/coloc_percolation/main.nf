process COLOC_PERCOLATION {
  label "process_medium"

  input:
    tuple path(annData), path(module_for_cp)
  
  output:
    // This marker changes only on actual execution, not on a cache hit.
    tuple path("post_coloc_percolation_ex_*_modules.rds"), path("refinement_run_id.txt"), emit: post_cp_module_rds

  script:
  def args = task.ext.args ?: ''
    """
    s11_coloc_percolation.R \
        ${args} \
        --annData ${annData} \
        --pph4 ${params.pph4_threshold} \
        --partly_connected_module ${module_for_cp} \
        --method ${params.module_refining_method}

    # Write the marker after successful refinement.
    python -c 'import uuid; print(uuid.uuid4())' > refinement_run_id.txt
    """
}