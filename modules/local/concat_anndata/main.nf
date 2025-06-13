process CONCAT_ANNDATA {
  tag "concat_anndata"
  label "process_medium"

  publishDir "${params.outdir}/results/anndata/", mode: params.publish_dir_mode, pattern:"*.h5ad"

  input:
    path(all_h5ad)
  
  output:
    path "coloc_input_anndata.h5ad", emit: full_anndata

  script:
  def args = task.ext.args ?: ''
    """
    export RETICULATE_PYTHON=\$(which python)
    
    ls *.h5ad > all_h5ad_input_list.txt
    
    s08_concat_anndata.R \
        ${args} \
        --input all_h5ad_input_list.txt \
        --output_file coloc_input_anndata.h5ad
    """

  stub:
    """
    touch complete_anndata.h5ad

    """
}
