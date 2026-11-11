process MASTER_TABLE {
  label "process_medium"

  input:
    path(annData)
    path(modules_fully_connected)
    // An empty list is valid when refinement produced no results.
    tuple path(refinement_rds, stageAs: 'refinement????/*'), val(refinement_run_id)

  output:
    tuple val(refinement_run_id), path("${params.coloc_id}_master_table_chunk_*.tsv"), emit: master_table_chunk

  script:
  def args = task.ext.args ?: ''
  def result_files = refinement_rds instanceof List ? refinement_rds : [refinement_rds]
  def result_list = result_files.collect { it.toString() }.join('\n')
  """
    # Build the manifest from staged inputs. With no refinement results, it contains only a blank line, which the R reader ignores.
    cat > post_cp_rds.txt <<'RDS_LIST'
${result_list}
RDS_LIST

    s12_create_master_table.R \
        ${args} \
        --annData ${annData} \
        --fully_connected_modules ${modules_fully_connected} \
        --post_cp_list post_cp_rds.txt \
        --chunk_size ${params.chunk_size} \
        --coloc_id ${params.coloc_id}
    """
}