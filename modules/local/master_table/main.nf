process MASTER_TABLE {
  label "process_medium"

  input:
    path(annData)
    path(modules_fully_connected)
    // An empty list is valid when refinement produced no results.
    path(refinement_rds, stageAs: 'refinement/*')

  output:
    path "${params.coloc_id}_master_table_chunk_*.tsv", emit: master_table_chunk

  script:
  def args = task.ext.args ?: ''
  """
    # Ensures the folder exists even when no refinement files were produced.
    mkdir -p refinement

    s12_create_master_table.R \
        ${args} \
        --annData ${annData} \
        --fully_connected_modules ${modules_fully_connected} \
        --post_cp_dir refinement \
        --chunk_size ${params.chunk_size} \
        --coloc_id ${params.coloc_id}
    """
}