process FIND_CS_OVERLAP_BY_CHR {
  tag "${params.coloc_id}_cs_overlap"
  label "process_high"

  input:
    tuple val(meta_chr_cs), path(all_cond_datasets_cs)
  
  output:
    tuple val(meta_chr_cs), path("${params.coloc_id}_chr${meta_chr_cs.chr_cs}_coloc_pairwise_guide_table.tsv"), optional:true, emit:coloc_pairwise_guide_table

  script:
  def args = task.ext.args ?: ''
    """
    s05_find_overlapping_cs.R \
        ${args} \
        --pipeline_path ${projectDir}/bin/ \
        --coloc_info_table ${all_cond_datasets_cs} \
        --chr_cs ${meta_chr_cs.chr_cs} \
        --coloc_id ${params.coloc_id}
    """
  
  stub:
    """
    touch ${params.coloc_id}_chr${meta_chr_cs.chr_cs}_coloc_pairwise_guide_table.tsv
    """
}