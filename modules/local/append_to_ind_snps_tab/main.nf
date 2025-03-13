#!/usr/bin/env nextflow

process APPEND_TO_IND_SNPS_TAB {
  
  label "process_single"
  
  // Define input
  input:
    tuple val(meta_study_id), path(ind_snps_table)
  
  // Define output
  output:
    path "${meta_study_id.study_id}_final_ind_snps_table.tsv"

  // Publish output file to specified directory   
  publishDir "results/gwas_and_loci_tables", mode:"copy", pattern:"*_final_ind_snps_table.tsv"

// Define the shell script to execute
  shell:
    '''
    cat !{ind_snps_table} >> !{meta_study_id.study_id}.tmp
    sort -r !{meta_study_id.study_id}.tmp | awk 'NR == 1 || $0 != prev { print; prev = $0 }' > !{meta_study_id.study_id}_final_ind_snps_table.tsv
    
    '''
}