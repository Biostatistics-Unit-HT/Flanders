process MAKE_COLOC_GUIDE_TABLE {
  tag "make_guide_table"
  label "process_medium"

  publishDir "${params.outdir}/results/anndata/", mode: params.publish_dir_mode, pattern:"*.csv"

  input:
    path(h5ad_file)
    path(previous_h5ad_studies)
    path(exclude_studies_file)
  
  output:
    path "coloc_guide_table.csv", emit: coloc_guide_table

  script:
  def args = task.ext.args ?: ''
  def exclude_oneside_opt = exclude_studies_file.name != 'NO_EXCLUDE_STUDIES' ? "--exclude_oneside ${exclude_studies_file}" : ''
  def exclude_bothsides_opt = previous_h5ad_studies.name != 'NO_PREVIOUS_H5AD_STUDIES' ? "--exclude_bothsides ${previous_h5ad_studies}" : ''
  
    """
    export RETICULATE_PYTHON=\$(which python)
    
    s09_make_coloc_guide.R \
        ${args} \
        --input ${h5ad_file} \
        ${exclude_bothsides_opt} \
        ${exclude_oneside_opt} \
        --output_file coloc_guide_table.csv
    """

  stub:
    """
    touch coloc_guide_table.csv

    """
}
