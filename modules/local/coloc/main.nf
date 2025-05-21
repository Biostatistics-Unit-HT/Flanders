process COLOC {
  tag "${params.coloc_id}"
  label "process_medium"
  
  input:
    tuple val(meta_chr_cs), path(coloc_pairs_by_batches), path(rds_files) // why not path(coloc_pairs_by_batches) ?
  
  output:
    path "${params.coloc_id}_chr${meta_chr_cs.chr_cs}_colocalization.table.all.tsv", emit:colocalization_table_all_by_chunk

  script:
  def args = task.ext.args ?: ''
    """
    s06_coloc.R \
        ${args} \
        --pipeline_path ${projectDir}/bin/ \
        --coloc_guide_table ${coloc_pairs_by_batches} \
        --chr_cs ${meta_chr_cs.chr_cs} \
        --coloc_id ${params.coloc_id}
    """

  stub:
    """
    touch ${params.coloc_id}_chr${meta_chr_cs.chr_cs}_colocalization.table.all.tsv
    """
}