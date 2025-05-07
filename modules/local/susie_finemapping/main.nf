process SUSIE_FINEMAPPING {
  tag "${meta_study_id.study_id}"
  label "process_high"
  
  publishDir "${params.outdir}/results/finemap/", mode: params.publish_dir_mode, pattern:"*_susie_finemap.rds", enabled: params.publish_susie

  input:
    tuple val(meta_study_id), val(meta_finemapping), path(bfile_dataset), val(meta_loci), path(gwas_final), path(gwas_final_index)
    val outdir
  
  output:
    path "*_susie_finemap.rds", optional:true, emit: susie_results_rds
    tuple val(meta_study_id), path ("*_cs_info_table.tsv"), optional:true, emit:susie_info_coloc_table
    tuple val(meta_study_id), val(meta_finemapping), val(meta_loci), path(gwas_final), path(gwas_final_index), path(bfile_dataset), path("failed_susie.txt"), optional:true, emit:failed_susie_loci
    path "*_switched_to_L1_finemap_report.tsv", optional:true, emit: switched_to_L1_finemap_report
    path "*_failed_finemap_report.tsv", optional:true, emit: failed_finemap_report

  script:
  def args = task.ext.args ?: ''
    """
    s04_susie_finemapping.R \
        ${args} \
        --pipeline_path ${projectDir}/bin/ \
        --chr ${meta_loci.chr} \
        --start ${meta_loci.start} \
        --end ${meta_loci.end} \
        --phenotype_id ${meta_loci.phenotype_id} \
        --dataset_aligned ${gwas_final} \
        --maf ${meta_finemapping.maf} \
        --bfile ${bfile_dataset[0].baseName} \
        --skip_dentist ${meta_finemapping.skip_dentist} \
        --cs_thresh ${meta_finemapping.cs_thresh} \
        --susie_max_iter ${params.susie_max_iters}\
        --results_path ${outdir} \
        --study_id ${meta_study_id.study_id}
    """

  stub:
    """
    touch ${meta_study_id.study_id}_${meta_loci.phenotype_id}_locus_chr${meta_loci.chr}_${meta_loci.start}_${meta_loci.end}_susie_finemap.rds
    touch ${meta_study_id.study_id}_${meta_loci.phenotype_id}_locus_chr${meta_loci.chr}_${meta_loci.start}_${meta_loci.end}_cs_info_table.tsv
    """
}