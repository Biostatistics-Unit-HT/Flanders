process MASTER_TABLE {
  conda '/ssu/bsssu/conda_envs_bsssu/flanders_r_2.3.11_patched'
  label "process_medium"

  input:
    path(annData)
    path(modules_fully_connected)
    path(rds_list_file_all)

  output:
    path "${params.coloc_id}_master_table_chunk_*.tsv", emit: master_table_chunk

  script:
  def args = task.ext.args ?: ''
  """
    s12_create_master_table.R \
        ${args} \
        --annData ${annData} \
        --fully_connected_modules ${modules_fully_connected} \
        --post_cp_list ${rds_list_file_all} \
        --chunk_size ${params.chunk_size} \
        --coloc_id ${params.coloc_id}
    """
}