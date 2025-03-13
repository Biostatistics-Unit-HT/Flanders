#!/usr/bin/env nextflow

process COJO_AND_FINEMAPPING {
  label "process_high"

  // Define input
  input:
    tuple val(meta_study_id), val(meta_finemapping), val(meta_loci), path(gwas_final), path(gwas_final_index), path(susie_error_message)
    val lauDir
  
  // Define output
  output:
    path "*_conditioned_loci.pdf", optional:true
    path "*_cojo_finemap.rds", optional:true
    tuple val(meta_study_id), path ("*_coloc_info_table.tsv"), optional:true, emit:cojo_info_coloc_table
    tuple val(meta_study_id), path ("*_ind_snps.tsv"), optional:true, emit:ind_snps_table

  // Publish output file to specified directory   
  publishDir "plots/", mode:"copy", pattern:"*_conditioned_loci.pdf"
  publishDir "results/finemap/", mode:"copy", pattern:"*_cojo_finemap.rds"
  
  // Tag the process with the study ID    
  tag "${meta_study_id.study_id}_cojo_finemap"

// Define the shell script to execute
  shell:
    '''
    Rscript --vanilla !{projectDir}/bin/s04_cojo_finemapping.R \
        --pipeline_path !{projectDir}/bin/ \
        --chr !{meta_loci.chr} \
        --start !{meta_loci.start} \
        --end !{meta_loci.end} \
        --phenotype_id !{meta_loci.phenotype_id} \
        --dataset_aligned !{gwas_final} \
        --p_thresh3 !{meta_finemapping.p_thresh3} \
        --maf !{meta_finemapping.maf} \
        --bfile !{meta_finemapping.bfile} \
        --skip_dentist !{meta_finemapping.skip_dentist} \
        --p_thresh4 !{meta_finemapping.p_thresh4} \
        --hole !{meta_finemapping.hole} \
        --cs_thresh !{meta_finemapping.cs_thresh} \
        --results_path !{lauDir} \
        --study_id !{meta_study_id.study_id}
    '''
}