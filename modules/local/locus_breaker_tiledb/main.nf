#!/usr/bin/env nextflow

process LOCUS_BREAKER_TILEDB {
  label "process_medium"
  conda '/ssu/gassu/conda_envs/scqtl'
  publishDir "${params.outdir}/results/gwas_and_loci_tables/", mode: params.publish_dir_mode


// Define input
  input:
  path(traits_list_table)

// Define output
  output:
    tuple path("*_interval.csv"), emit:locus_breaker_tdb_intervals
    tuple path("dummy_index"), path("*_segment.csv"), emit:locus_breaker_tdb_segments

// Define the shell script to execute
  script:
    """
    scqtl --workers 4 \
      export \
      --table ${traits_list_table} \
      --uri-path ${params.tiledb_uri} \
      --out_lb ${params.path_tiledb_lb_out} \
      --maf ${params.tiledb_lb_maf} \
      --type-sumstat ${params.tiledb_lb_typesumstat} \
      --pvalue-sig ${params.tiledb_lb_pvalue_sig} \
      --pvalue-limit ${params.tiledb_lb_pvalue_limit} \
      --hole ${params.tiledb_lb_hole} \
      --locus-max-size ${params.tiledb_large_locus_size} \
      --locusbreaker
    
    touch dummy_index
    """
}

