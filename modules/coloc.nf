#!/usr/bin/env nextflow

process COLOC {

  // Define input
  input:
    tuple val(meta_chr_cs), path(coloc_pairs_by_batches) // why not path(coloc_pairs_by_batches) ?
  
  // Define output
  output:
    path "${params.coloc_id}_chr${meta_chr_cs.chr_cs}_colocalization.table.all.tsv", emit:colocalization_table_all_by_chunk

  // Tag the process with the study ID    
  tag "${params.coloc_id}_coloc"

// Define the shell script to execute
  shell:
    '''
    Rscript --vanilla !{projectDir}/bin/s06_coloc.R \
        --pipeline_path !{projectDir}/bin/ \
        --coloc_guide_table !{coloc_pairs_by_batches} \
        --chr_cs !{meta_chr_cs.chr_cs} \
        --coloc_id !{params.coloc_id}
    '''
}