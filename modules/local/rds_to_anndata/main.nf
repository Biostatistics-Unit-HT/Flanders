process RDS_TO_ANNDATA {
  tag "rds_to_anndata"
  label "process_high"

  publishDir "${params.outdir}/results/anndata/", mode: params.publish_dir_mode, pattern:"*.h5ad"

  input:
    path(all_rds)
  
  output:
    path "*.h5ad", emit: finemap_anndata

  script:
  def args = task.ext.args ?: ''
    """
    export RETICULATE_PYTHON=\$(which python)
    
    ls *.rds > all_rds_input_list.txt
    
    s07_anndata_concat.R \
        ${args} \
        --input all_rds_input_list.txt \
        --output_file finemap_results.h5ad
    """

  stub:
    """
    touch finemap_results.h5ad

    """
}
