process MAKE_COLOC_GUIDE_TABLE {
  tag "make_guide_table"
  label "process_medium"

  publishDir "${params.outdir}/results/anndata/", mode: params.publish_dir_mode, pattern:"*.csv"

  input:
    path(h5ad_file)
    path(studies_to_exclude)
  
  output:
    path "coloc_guide_table.csv", emit: coloc_guide_table

  script:
  def args = task.ext.args ?: ''
    """
    export RETICULATE_PYTHON=\$(which python)
    
    s09_make_coloc_guide.R \
        ${args} \
        --input ${h5ad_file} \
        --exclude ${studies_to_exclude} \
        --output_file coloc_guide_table.csv
    """

  stub:
    """
    touch coloc_guide_table.csv

    """
}
