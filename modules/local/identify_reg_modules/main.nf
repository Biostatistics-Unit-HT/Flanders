process IDENTIFY_REG_MODULES {
  tag "${params.coloc_id}_identify_reg_modules"
  label "process_low"
  cache 'deep'
  
  input:
    path(coloc_results_all)

  output:
    path "modules_fully_connected.tsv", optional:true, emit: modules_fully_connected
    path "partly_connected_modules_*.tsv", optional:true, emit: modules_partly_connected

  script:
  def args = task.ext.args ?: ''
    """
    s10_identify_reg_modules.R \
        ${args} \
        --coloc_results ${coloc_results_all} \
        --pph4 ${params.pph4_threshold} \
        --coloc_id ${params.coloc_id}
    """
}