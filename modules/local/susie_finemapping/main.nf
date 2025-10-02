process SUSIE_FINEMAPPING {
  tag "${meta_study_id.study_id}"
  label "process_high"
  
  publishDir "${params.outdir}/results/finemap/", mode: params.publish_dir_mode, pattern:"*_anndata.h5ad"

  input:
    // make GWAS/index inputs optional so the process can be invoked with placeholders
    tuple val(meta_study_id), val(meta_finemapping), path(bfile_dataset), val(meta_loci), path(gwas_final), path(gwas_final_index)
    val outdir

  output:
    tuple val(meta_study_id), path ("*_cs_info_table.tsv"), optional:true, emit:susie_info_coloc_table
    path "*_anndata.h5ad", optional:true, emit: susie_results_h5ad
    path "*_FINEMAPPED_L1_prior_variance_too_large.tsv", optional:true, emit: finemapped_L1_prior_variance_too_large
    path "*_FINEMAPPED_L1_IBSS_algorithm_did_not_converge.tsv", optional:true, emit: finemapped_L1_IBSS_algorithm_did_not_converge
    path "*_FINEMAPPED_L1_recover_after_susie_QC.tsv", optional:true, emit: finemapped_L1_recover_after_susie_QC
    path "*_NOT_FINEMAPPED_no_credible_sets_found.tsv", optional:true, emit: not_finemapped_no_credible_sets_found
    path "*_NOT_FINEMAPPED_no_variants_from_locus_in_LD_ref.tsv", optional:true, emit: not_finemapped_no_variants_from_locus_in_LD_ref

  script:
    /*
      Defensive binding:
      - When params.tiledb is true we won't pass an index flag to the R script.
      - When gwas_final / gwas_final_index are not provided (null), create/use safe placeholders
        so Nextflow does not complain about null path values and the R script receives valid args.
    */
    def args = task.ext.args ?: ''
    // Robustly locate interval and segment files whether gwas_final is a List or single Path,
    // and whether the segment might have been bound to gwas_final_index.
    def gwas_segment = null
    def gwas_interval = null
    if (params.tiledb) {
      if (gwas_final instanceof List) {
        gwas_segment = gwas_final.find { it.name.contains('segment') }
        gwas_interval = gwas_final.find { it.name.contains('interval') }
      } else {
        // If gwas_final is a single Path, interval may be gwas_final and segment may be gwas_final_index
        gwas_interval = gwas_final
        gwas_segment = gwas_final_index
      }
    } else {
      // non-tiledb: dataset_aligned is gwas_final
      gwas_segment = gwas_final
    }

    // Error if segment not found
    if (params.tiledb && !gwas_segment) {
      // print helpful context and fail fast rather than passing "null" into R
      throw new IllegalStateException("SUSIE_FINEMAPPING: couldn't find GWAS 'segment' file. gwas_final=${gwas_final}, gwas_final_index=${gwas_final_index}")
    }

    def batch_arg = params.tiledb && gwas_interval ? "--batch ${gwas_interval}" : ""
    


    // Build index param defensively
    def index_param = ""
    if (!params.tiledb) {
      if (gwas_final_index) {
        index_param = "--index_file ${gwas_final_index}"
      } else {
        // will create dummy_index file in the work dir at runtime
        index_param = "--index_file dummy_index"
      }
    }

    // In case bfile_dataset is a path tuple/list, keep the same access used previously
    def bfile_base = ''
    try {
      bfile_base = bfile_dataset[0].baseName
    } catch (all) {
      // fallback if bfile_dataset is a single path
      bfile_base = bfile_dataset.baseName
    }

    """
    # If necessary, create a dummy index file so the index path is non-null on-disk
    if [ "${index_param}" = "--index_file dummy_index" ]; then
      touch dummy_index
    fi

    export RETICULATE_PYTHON=\$(which python)

    s04_susie_finemapping.R \
        ${args} \
        --pipeline_path ${projectDir}/bin/ \
        --chr ${meta_loci.chr} \
        --start ${meta_loci.start} \
        --end ${meta_loci.end} \
        --phenotype_id ${meta_loci.phenotype_id} \
        --dataset_aligned ${gwas_segment} \
        --maf ${meta_finemapping.maf} \
        --bfile ${bfile_dataset[0].baseName} \
        --skip_dentist ${meta_finemapping.skip_dentist} \
        --cs_thresh ${meta_finemapping.cs_thresh} \
        --susie_max_iter ${params.susie_max_iter}\
        --susie_qc_cs_bf_thr ${params.susie_qc_cs_bf_thr}\
        --susie_qc_pval_thr ${params.susie_qc_pval_thr}\
        --susie_qc_mean_r2_thr ${params.susie_qc_mean_r2_thr}\
        --susie_qc_min_r2_thr ${params.susie_qc_min_r2_thr}\
        --publish_susie ${params.publish_susie}\
        --results_path ${outdir} \
        --study_id ${meta_study_id.study_id} \
        ${batch_arg}
    """
  stub:
    """
    touch ${meta_study_id.study_id}_${meta_loci.phenotype_id}_locus_chr${meta_loci.chr}_${meta_loci.start}_${meta_loci.end}_cs_info_table.tsv
    touch batch${gwas_interval}_anndata.h5ad
    """
}