process APPEND_TO_IND_SNPS_TAB {
  tag "${meta_study_id.study_id}"
  label "process_single"

  publishDir "${params.outdir}/results/gwas_and_loci_tables", mode: params.publish_dir_mode, pattern:"*_final_ind_snps_table.tsv"
  
  input:
    tuple val(meta_study_id), path(ind_snps_table)
  
  output:
    path "${meta_study_id.study_id}_final_ind_snps_table.tsv"  

  script:
    """
    cat ${ind_snps_table} >> ${meta_study_id.study_id}.tmp
    sort -r ${meta_study_id.study_id}.tmp | awk 'NR == 1 || \$0 != prev { print; prev = \$0 }' > ${meta_study_id.study_id}_final_ind_snps_table.tsv
    """

  stub:
    """
    touch ${meta_study_id.study_id}_final_ind_snps_table.tsv
    """
}