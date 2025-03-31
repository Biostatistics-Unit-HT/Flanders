#!/usr/bin/env nextflow
process MUNG_AND_LOCUS_BREAKER {
  tag "${meta_study_id.study_id}_mung_and_locus_break"
  label "process_high"

  // TODO: This works but it's slow for the moment we use a global conda env for the whole pipeline
  // conda "${moduleDir}/environment.yml"

  // Publish output file to specified directory   
  publishDir "${params.outdir}/results/gwas_and_loci_tables", mode: params.publish_dir_mode, pattern:"${meta_study_id.study_id}_dataset_aligned.tsv.gz"
  publishDir "${params.outdir}/results/gwas_and_loci_tables", mode: params.publish_dir_mode, pattern:"${meta_study_id.study_id}_loci.tsv"  
  
  // Define input - tuple (similar to a list) of: study id, other specific metadata parameters, gwas sum stats
  input:
    tuple val(meta_study_id), val(meta_parameters), path(gwas_input)
    val (validation)
    path chain_file

  // Define output - keep carrying the study id (will be used for all downstream processes!), output munged gwas .rds
  output:
    tuple val(meta_study_id), path("${meta_study_id.study_id}_dataset_aligned_indexed.tsv.gz"), path("${meta_study_id.study_id}_dataset_aligned_indexed.tsv.gz.tbi"), emit:dataset_munged_aligned
    path("${meta_study_id.study_id}_loci_NO_HLA.tsv"), optional:true, emit:loci_table
    path("${meta_study_id.study_id}_dataset_aligned.tsv.gz")
    path("${meta_study_id.study_id}_loci.tsv"), optional:true

  script:
    """
    s02_sumstat_munging_and_aligning.R \
        --input ${gwas_input} \
        --is_molQTL ${meta_parameters.is_molQTL} \
        --run_liftover ${meta_parameters.run_liftover} \
        --key ${meta_parameters.key} \
        --rsid_lab ${meta_parameters.rsid_lab} \
        --chr_lab ${meta_parameters.chr_lab} \
        --pos_lab ${meta_parameters.pos_lab} \
        --a1_lab ${meta_parameters.a1_lab} \
        --a0_lab ${meta_parameters.a0_lab} \
        --effect_lab ${meta_parameters.effect_lab} \
        --se_lab ${meta_parameters.se_lab} \
        --freq_lab ${meta_parameters.freq_lab} \
        --pvalue_lab ${meta_parameters.pvalue_lab} \
        --n_lab ${meta_parameters.n_lab} \
        --type ${meta_parameters.type} \
        --sdY ${meta_parameters.sdY} \
        --s ${meta_parameters.s} \
        --bfile ${meta_parameters.bfile} \
        --grch ${meta_parameters.grch} \
        --maf ${meta_parameters.maf} \
        --p_thresh1 ${meta_parameters.p_thresh1} \
        --p_thresh2 ${meta_parameters.p_thresh2} \
        --hole ${meta_parameters.hole} \
        --study_id ${meta_study_id.study_id} \
        --threads ${task.cpus}
    """
}