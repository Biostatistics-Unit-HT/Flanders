process GET_STUDY_FROM_ANNDATA {
  tag "get_study_from_anndata"
  label "process_medium"

  input:
    path(all_h5ad)
  
  output:
    path "previous_h5ad_studies.tsv", emit: previous_study_ids

  script:
  def args = task.ext.args ?: ''
    """
    export RETICULATE_PYTHON=\$(which python)
    
    ls *.h5ad > all_h5ad_input_list.txt
    
    get_study_and_pheno_id.py \
        ${args} \
        --input all_h5ad_input_list.txt \
        --output_file previous_h5ad_studies.tsv
    """

  stub:
    """
    touch complete_anndata.h5ad

    """
}
