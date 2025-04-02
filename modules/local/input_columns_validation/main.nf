process INPUT_COLUMNS_VALIDATION {
    label "process_single"

    input:
        path(file_list)
        val(base_dir)

    output:
        path file_list, emit: table_out

    script:
    def args = task.ext.args ?: ''
    def liftover_flag = params.run_liftover ? "--liftover" : ""
    """
    validate_columns.py --launchdir ${base_dir} --table ${file_list} $liftover_flag $args
    
    if [[ \$? -ne 0 ]]; then
        echo "ERROR: Validation failed. Check logs for details." >&2
        exit 1
    fi
    """

    stub:
    def liftover_flag = params.run_liftover ? "--liftover" : ""
    """
    echo "${liftover_flag}" > ok.txt
    """
}
