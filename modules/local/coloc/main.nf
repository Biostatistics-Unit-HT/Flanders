process COLOC {
  tag "${params.coloc_id}"
  label "process_medium"
  
  // Define input
  input:
    tuple file(annData), val(coloc_pairs_by_batches)
  
  // Define output
  output:
    path "${params.coloc_id}_colocalization.table.all.tsv", emit:colocalization_chunk

  // Tag the process with the study ID    
  tag "${params.coloc_id}_coloc"

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