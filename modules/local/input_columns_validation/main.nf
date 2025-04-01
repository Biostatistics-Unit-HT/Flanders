process INPUT_COLUMNS_VALIDATION {
    label "process_single"

    input:
        val(file_list)

    output:
        val(file_list), emit: table_out

    script:
    def liftover_flag = params.run_liftover ? "--liftover" : ""
    """
    python ${projectDir}/bin/validate_columns.py --launchdir ${launchDir} --table ${file_list} $liftover_flag
    
    if [[ \$? -ne 0 ]]; then
        echo "ERROR: Validation failed. Check logs for details." >&2
        exit 1
    fi

    """
}
