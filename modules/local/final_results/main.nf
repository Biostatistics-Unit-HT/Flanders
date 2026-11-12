process FINAL_RESULTS {
  label "process_low"

  publishDir "${params.outdir}/results/master_table", mode: "copy", pattern:"${params.coloc_id}_modules_JF.h5ad"
  publishDir "${params.outdir}/results/master_table", mode: "copy", pattern:"${params.coloc_id}_master_table.tsv"

  input:
    tuple path(annData), path(joint_result_rds, stageAs: 'joint_results/*')

  output:
    path "${params.coloc_id}_modules_JF.h5ad"
    path "${params.coloc_id}_master_table.tsv"

  script:
  def args = task.ext.args ?: ''
    """
    # Ensure the input directory exists.
    mkdir -p joint_results

    # Select Python from the task's container or Conda environment.
    export RETICULATE_PYTHON=\$(which python)

    s14_output_final_results.R \
        ${args} \
        --joint_finemapping_dir joint_results \
        --annData ${annData} \
        --coloc_id ${params.coloc_id}
    """
}