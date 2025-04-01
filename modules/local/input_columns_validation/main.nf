process INPUT_COLUMNS_VALIDATION {
    label "process_single"

    input:
        val(file_list)

    output:
        val(file_list), emit: table_out

    script:
    """
    python ${projectDir}/bin/validate_columns.py --launchdir ${launchDir} --liftover ${params.run_liftover} --table ${file_list}
    
    if [[ \$? -ne 0 ]]; then
        echo "ERROR: Validation failed. Check logs for details." >&2
        exit 1
    fi

    """
}
