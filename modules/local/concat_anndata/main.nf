process CONCAT_ANNDATA {
    tag "concat_anndata"
    label "process_medium"

    publishDir "${params.outdir}/results/anndata/", mode: params.publish_dir_mode, pattern:"*.h5ad"

    input:
    path(all_h5ad, stageAs: 'input_file???/*')
    val(output_name)
    
    output:
    path "${output_name}", emit: anndata

    script:
    def args = task.ext.args ?: ''

    """
    export RETICULATE_PYTHON=\$(which python)
    
    ls input_file*/*.h5ad > all_h5ad_input_list.txt
    
    s08_concat_anndata.R \
        ${args} \
        --input all_h5ad_input_list.txt \
        --output_file ${output_name}
    """

    stub:
    """
    touch ${output_name}
    """
}