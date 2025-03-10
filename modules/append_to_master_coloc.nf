#!/usr/bin/env nextflow

process APPEND_TO_MASTER_COLOC {
  
  label 'appending_tables'

  // Define input
  input:
    tuple val(meta_study_id), path(info_coloc_table)
  
  // Define output
  output:
    path "${meta_study_id.study_id}_coloc_info_master_table.tsv"

  // Publish output file to specified directory   
  publishDir "results/coloc_info_tables", mode:"copy", pattern:"*_coloc_info_master_table.tsv"

// Define the shell script to execute
  shell:
    '''
    echo -e "study_id\tphenotype_id\tcredible_set\ttop_pvalue\tpath_rds\tpath_ind_snps\tchr" > !{meta_study_id.study_id}_coloc_info_master_table.tsv
    cat !{info_coloc_table} >> !{meta_study_id.study_id}_coloc_info_master_table.tsv
    
    '''
}