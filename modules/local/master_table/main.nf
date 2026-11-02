process MASTER_TABLE {
  conda '/ssu/bsssu/conda_envs_bsssu/flanders_r_2.3.8'
  label "process_medium"

  publishDir "${params.outdir}/results/master_table", mode: "copy", pattern:"${params.coloc_id}_master_table.tsv"

  input:
    path(annData)
    path(modules_fully_connected)
    path(rds_list_file_all)

  output:
    path "${params.coloc_id}_master_table_chunk_*.tsv", emit: master_table_chunk
    path "${params.coloc_id}_master_table.tsv", emit: master_table

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