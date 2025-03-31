#!/usr/bin/env nextflow

process MUNG_AND_ALIGN {

  // Define input - tuple (similar to a list) of: study id, other specific metadata parameters, gwas sum stats
  input:
    tuple val(meta_study_id), val(meta_parameters), path(gwas_input)

  // Define output - keep carrying the study id (will be used for all downstream processes!), output munged gwas .rds
  output:
    tuple val(meta_study_id), path("${meta_study_id.study_id}_dataset_aligned.tsv.gz"), emit:dataset_munged_aligned

  // Publish output file to specified directory   
  publishDir "${params.outdir}/results/gwas_and_loci_tables/", mode:"copy", pattern:"${meta_study_id.study_id}_dataset_aligned.tsv.gz"
   
  // Tag the process with the study ID    
  tag "${meta_study_id.study_id}_mung"

// Define the shell script to execute - IMPORTANT! Use 3 single apostrophes, otherwise bash variables will not be recognised
  shell:
    '''
    Rscript --vanilla !{projectDir}/bin/s02_sumstat_munging_and_aligning.R \
        --pipeline_path !{projectDir}/bin/ \
        --input !{gwas_input} \
        --is_molQTL !{meta_parameters.is_molQTL} \
        --key !{meta_parameters.key} \
        --rsid_lab !{meta_parameters.rsid_lab} \
        --chr_lab !{meta_parameters.chr_lab} \
        --pos_lab !{meta_parameters.pos_lab} \
        --a1_lab !{meta_parameters.a1_lab} \
        --a0_lab !{meta_parameters.a0_lab} \
        --effect_lab !{meta_parameters.effect_lab} \
        --se_lab !{meta_parameters.se_lab} \
        --freq_lab !{meta_parameters.freq_lab} \
        --pvalue_lab !{meta_parameters.pvalue_lab} \
        --n_lab !{meta_parameters.n_lab} \
        --type !{meta_parameters.type} \
        --sdY !{meta_parameters.sdY} \
        --s !{meta_parameters.s} \
        --bfile !{meta_parameters.bfile} \
        --grch !{meta_parameters.grch} \
        --study_id !{meta_study_id.study_id}  
    '''
}