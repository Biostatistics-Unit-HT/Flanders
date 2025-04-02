process LOCUS_BREAKER {
  tag "${meta_study_id.study_id}_locus_b"
  label "process_low"

  publishDir "${params.outdir}/results/gwas_and_loci_tables/", mode: params.publish_dir_mode, pattern:"${meta_study_id.study_id}_loci.tsv"

  input:
    tuple val(meta_study_id), val(meta_locus_breaker), path(dataset_munged_aligned)
  
  output:
    path "${meta_study_id.study_id}_loci.tsv", optional:true, emit:loci_table

  script:
  def args = task.ext.args ?: ''
    """
    s03_locus_breaker.R \
        ${args} \
        --pipeline_path ${projectDir}/bin/ \
        --dataset_aligned ${dataset_munged_aligned} \
        --maf ${meta_locus_breaker.maf} \
        --p_thresh1 ${meta_locus_breaker.p_thresh1} \
        --p_thresh2 ${meta_locus_breaker.p_thresh2} \
        --hole ${meta_locus_breaker.hole} \
        --study_id ${meta_study_id.study_id}
    """

  stub:
    """
    touch ${meta_study_id.study_id}_loci.tsv
    """
}