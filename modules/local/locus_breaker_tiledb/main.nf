#!/usr/bin/env nextflow

process LOCUS_BREAKER_TILEDB {
  label "process_multi"
  //conda '/ssu/gassu/conda_envs/tdbsumstat'
  conda '/software/cardinal_analysis/ht/conda_envs/tdbsumstat'

  publishDir "${params.outdir}/results/gwas_and_loci_tables/", mode: params.publish_dir_mode


// Define input
  input:
  tuple  val(batch_index), path(traits_list_table)

// Define output
  output:
    path("*_interval.csv"), emit:locus_breaker_tdb_intervals, optional: true
    tuple path("dummy_index"), path("*_segment.csv"), emit:locus_breaker_tdb_segments, optional: true

// Define the shell script to execute
  script:
    """
    tdbsumstat \
      export \
      --table-lb ${traits_list_table} \
      --uri-path ${params.tiledb_uri} \
      --out ${params.path_tiledb_lb_out} \
      --maf-lb ${params.tiledb_lb_maf} \
      --type-sumstat ${params.tiledb_lb_typesumstat} \
      --hole-lb ${params.tiledb_lb_hole} \
      --locus-max-size-lb ${params.tiledb_large_locus_size} \
      --locusbreaker \
      --batch-name ${batch_index}
    
    touch dummy_index
    """
}

