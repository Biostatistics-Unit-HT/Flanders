process JOINT_FINEMAP_BY_CHUNK {
  label "process_medium"

  input:
    tuple path(annData), path(master_table_by_chunk), val(refinement_run_id)

  output:
    path "${params.coloc_id}_joint_finemapping_master_table_chunk_*.tsv", emit: joint_finemapping_tsv_chunk
    tuple val(refinement_run_id), path("${params.coloc_id}_joint_finemapping_lABF_matrix_chunk_*.rds"), emit: joint_finemapping_rds_chunk

  script:
  def args = task.ext.args ?: ''
  """
    s13_joint_finemapping.R \
        ${args} \
        --annData ${annData} \
        --master_table_chunk ${master_table_by_chunk} \
        --coloc_id ${params.coloc_id}
    """
}