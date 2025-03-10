#!/usr/bin/env nextflow

process SUSIE_FINEMAPPING {

  // Define input
  input:
    tuple val(meta_study_id), val(meta_finemapping), path(gwas_final)
    val lauDir
  
  // Define output
  output:
    path "*_susie_finemap.rds", optional:true
    tuple val(meta_study_id), path ("*_coloc_info_table.tsv"), optional:true, emit:susie_info_coloc_table
    tuple val(meta_study_id), path ("*_ind_snps.tsv"), optional:true, emit:ind_snps_table    
    tuple val(meta_study_id), val(meta_finemapping), path(gwas_final), path("failed_susie.txt"), optional:true, emit:failed_susie_loci

  // Publish output file to specified directory   
  publishDir "results/finemap/", mode:"copy", pattern:"*_susie_finemap.rds"

  // Tag the process with the study ID    
  tag "${meta_study_id.study_id}_susie_finemap"

// Define the shell script to execute
  shell:
    '''
    Rscript --vanilla !{projectDir}/bin/s04_susie_finemapping.R \
        --pipeline_path !{projectDir}/bin/ \
        --chr !{meta_finemapping.chr} \
        --start !{meta_finemapping.start} \
        --end !{meta_finemapping.end} \
        --phenotype_id !{meta_finemapping.phenotype_id} \
        --dataset_aligned !{gwas_final} \
        --maf !{meta_finemapping.maf} \
        --bfile !{meta_finemapping.bfile} \
        --skip_dentist !{meta_finemapping.skip_dentist} \
        --cs_thresh !{meta_finemapping.cs_thresh} \
        --results_path !{lauDir} \
        --study_id !{meta_study_id.study_id}
    '''
}