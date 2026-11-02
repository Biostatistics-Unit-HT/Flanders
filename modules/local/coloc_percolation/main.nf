process COLOC_PERCOLATION {
  conda '/ssu/bsssu/conda_envs_bsssu/flanders_r_2.3.8'
  label "process_medium"

  input:
    tuple path(annData), path(module_for_cp)
  
  output:
    path "post_coloc_percolation_ex_*_modules.rds", emit: post_cp_module_rds

  script:
  def args = task.ext.args ?: ''
    """
    s11_coloc_percolation.R \
        ${args} \
        --annData ${annData} \
        --pph4 ${params.pph4_threshold} \
        --partly_connected_module ${module_for_cp} \
        --method ${params.module_refining_method}
    """
}