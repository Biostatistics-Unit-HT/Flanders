process COLOC {
  tag "${coloc_id}"
  label "process_medium"
  
  input:
    tuple path(annData), path(coloc_pairs_by_batches)
    val(coloc_id)
  
  output:
    path "${coloc_id}_colocalization.table.all.tsv", emit:colocalization_chunk

  script:
  def args = task.ext.args ?: ''
    """
    s06_coloc.R \
      ${args} \
      --annData ${annData} \
      --coloc_guide_table ${coloc_pairs_by_batches} \
      --coloc_id ${coloc_id}
    """

  stub:
    """
    touch ${coloc_id}_colocalization.table.all.tsv
    """
}