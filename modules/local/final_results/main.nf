process FINAL_RESULTS {
  label "process_low"

  publishDir "${params.outdir}/results/master_table", mode: "copy", pattern:"${params.coloc_id}_modules_JF.h5ad"
  publishDir "${params.outdir}/results/master_table", mode: "copy", pattern:"${params.coloc_id}_master_table.tsv"

  input:
    tuple path(annData), path(joint_result_rds, stageAs: 'joint????/*'), val(refinement_run_id)

  output:
    path "${params.coloc_id}_modules_JF.h5ad"
    path "${params.coloc_id}_master_table.tsv"

  script:
  def args = task.ext.args ?: ''
  def result_files = joint_result_rds instanceof List ? joint_result_rds : [joint_result_rds]
  def result_list = result_files.collect { it.toString() }.join('\n')
    """
    # List the staged RDS inputs, so both staging and caching track them.
    cat > all_joint_finemapping_rds.txt <<'RDS_LIST'
${result_list}
RDS_LIST

    # Select Python from the task's container or Conda environment.
    export RETICULATE_PYTHON=\$(which python)
    
    s14_output_final_results.R \
        ${args} \
        --all_joint_finemapping_obj all_joint_finemapping_rds.txt \
        --annData ${annData} \
        --coloc_id ${params.coloc_id}
    """
}