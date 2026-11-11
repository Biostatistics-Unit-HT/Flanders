process FINAL_RESULTS {
  label "process_low"

  publishDir "${params.outdir}/results/master_table", mode: "copy", pattern:"${params.coloc_id}_modules_JF.h5ad"
  publishDir "${params.outdir}/results/master_table", mode: "copy", pattern:"${params.coloc_id}_master_table.tsv"

  input:
    tuple path(annData), path(all_joint_finemapping_obj)

  output:
    path "${params.coloc_id}_modules_JF.h5ad"
    path "${params.coloc_id}_master_table.tsv"

  script:
  def args = task.ext.args ?: ''
    """
    # Select Python from the task's container or Conda environment.
    export RETICULATE_PYTHON=\$(which python)

    s14_output_final_results.R \
        ${args} \
        --all_joint_finemapping_obj ${all_joint_finemapping_obj} \
        --annData ${annData} \
        --coloc_id ${params.coloc_id}
    """
}