process INPUT_COLUMNS_VALIDATION {
    label "process_single"

    input:
        path(file_list)
        val(base_dir)

    output:
        path file_list, emit: table_out

    script:
    def liftover_flag = params.run_liftover ? "--liftover" : ""
    """
    validate_columns.py --launchdir ${base_dir} --table ${file_list} $liftover_flag
    
    if [[ \$? -ne 0 ]]; then
        echo "ERROR: Validation failed. Check logs for details." >&2
        exit 1
    fi

    """
}
