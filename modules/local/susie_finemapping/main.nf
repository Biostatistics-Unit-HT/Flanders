process SUSIE_FINEMAPPING {
  tag "${meta_study_id.study_id}_susie_finemap"
  label "process_high"
  
  publishDir "${params.outdir}/results/finemap/", mode: params.publish_dir_mode, pattern:"*_susie_finemap.rds"

  input:
    tuple val(meta_study_id), val(meta_finemapping), val(meta_loci), path(gwas_final), path(gwas_final_index)
    val lauDir
  
  output:
    path "*_susie_finemap.rds", optional:true
    tuple val(meta_study_id), path ("*_coloc_info_table.tsv"), optional:true, emit:susie_info_coloc_table
    tuple val(meta_study_id), path ("*_ind_snps.tsv"), optional:true, emit:ind_snps_table    
    tuple val(meta_study_id), val(meta_finemapping), val(meta_loci), path(gwas_final), path(gwas_final_index), path("failed_susie.txt"), optional:true, emit:failed_susie_loci

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
        --bfile ${meta_finemapping.bfile} \
        --skip_dentist ${meta_finemapping.skip_dentist} \
        --cs_thresh ${meta_finemapping.cs_thresh} \
        --results_path ${lauDir} \
        --study_id ${meta_study_id.study_id}
    """

  stub:
    """
    touch ${meta_study_id.study_id}_${meta_loci.phenotype_id}_locus_chr${meta_loci.chr}_${meta_loci.start}_${meta_loci.end}_susie_finemap.rds
    touch ${meta_study_id.study_id}_${meta_loci.phenotype_id}_locus_chr${meta_loci.chr}_${meta_loci.start}_${meta_loci.end}_susie_coloc_info_table.tsv
    touch ${meta_study_id.study_id}_${meta_loci.phenotype_id}_locus_chr${meta_loci.chr}_${meta_loci.start}_${meta_loci.end}_ind_snps.tsv
    """
}