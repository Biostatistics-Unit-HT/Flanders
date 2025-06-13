process GET_STUDY_FROM_ANNDATA {
  tag "get_study_from_anndata"
  label "process_medium"

  publishDir "${params.outdir}/results/anndata/", mode: params.publish_dir_mode, pattern:"*.tsv"

  input:
    path(all_h5ad)
  
  output:
    path "studies_excluded_previous_h5ad.tsv", emit: previous_study_ids

  script:
  def args = task.ext.args ?: ''
    """
    export RETICULATE_PYTHON=\$(which python)
    
    ls *.h5ad > all_h5ad_input_list.txt
    
    get_study_and_pheno_id.py \
        ${args} \
        --input all_h5ad_input_list.txt \
        --output_file studies_excluded_previous_h5ad.tsv
    """

  stub:
    """
    touch complete_anndata.h5ad

    """
}
