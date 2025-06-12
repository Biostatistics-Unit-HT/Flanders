process COLOC {
  tag "${params.coloc_id}"
  label "process_medium"
  
  input:
    tuple file(annData), val(coloc_pairs_by_batches)
  
  output:
    path "${params.coloc_id}_colocalization.table.all.tsv", emit:colocalization_chunk

  script:
  def args = task.ext.args ?: ''
    """
    s06_coloc.R \
      ${args} \
        --annData ${annData} \
        --coloc_guide_table ${coloc_pairs_by_batches} \
        --coloc_id ${params.coloc_id}
    """

  stub:
    """
    touch ${params.coloc_id}_colocalization.table.all.tsv
    """
}