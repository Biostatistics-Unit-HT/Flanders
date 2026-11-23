process COLOC_PERCOLATION {
  label "process_medium"

  input:
    tuple path(annData), path(module_for_cp)
    val(pph4_threshold)
    val(module_refining_method)
  
  output:
    path "post_coloc_percolation_ex_*_modules.rds", emit: post_cp_module_rds

  script:
  def args = task.ext.args ?: ''
    """
    s11_coloc_percolation.R \
        ${args} \
        --annData ${annData} \
        --pph4 ${pph4_threshold} \
        --partly_connected_module ${module_for_cp} \
        --method ${module_refining_method}
    """
}