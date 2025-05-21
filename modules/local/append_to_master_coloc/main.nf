process APPEND_TO_MASTER_COLOC {
  tag "${meta_study_id.study_id}"
  label "process_single"

  publishDir "${params.outdir}/results/coloc_info_tables", mode: params.publish_dir_mode, pattern:"*_coloc_info_master_table.tsv"

  input:
    tuple val(meta_study_id), path(info_coloc_table)
  
  output:
    path "${meta_study_id.study_id}_coloc_info_master_table.tsv", emit: coloc_master

  script:
    """
    echo -e "credible_set_name\tcredible_set_snps\tstudy_id\tphenotype_id\tchr\tstart\tend\ttop_pvalue\tpath_rds\tsnp\ta1\ta0\tfreq\tN\tbC\tbC_se" > ${meta_study_id.study_id}_coloc_info_master_table.tsv
    cat ${info_coloc_table} >> ${meta_study_id.study_id}_coloc_info_master_table.tsv
    """

  stub:
    """
    touch ${meta_study_id.study_id}_coloc_info_master_table.tsv
    """
}