#!/usr/bin/env nextflow

process LOCUS_BREAKER {
  label "process_low"

  // Define input
  input:
    tuple val(meta_study_id), val(meta_locus_breaker), path(dataset_munged_aligned)
  
  // Define output
  output:
    path "${meta_study_id.study_id}_loci.tsv", optional:true, emit:loci_table

  // Publish output file to specified directory   
  publishDir "${params.outdir}/results/gwas_and_loci_tables/", mode:"copy", pattern:"${meta_study_id.study_id}_loci.tsv"
  
  // Tag the process with the study ID    
  tag "${meta_study_id.study_id}_locus_b"

// Define the shell script to execute
  shell:
    '''
    Rscript --vanilla !{projectDir}/bin/s03_locus_breaker.R \
        --pipeline_path !{projectDir}/bin/ \
        --dataset_aligned !{dataset_munged_aligned} \
        --maf !{meta_locus_breaker.maf} \
        --p_thresh1 !{meta_locus_breaker.p_thresh1} \
        --p_thresh2 !{meta_locus_breaker.p_thresh2} \
        --hole !{meta_locus_breaker.hole} \
        --study_id !{meta_study_id.study_id}  
    '''
}