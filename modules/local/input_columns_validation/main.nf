process INPUT_COLUMNS_VALIDATION {
    label "process_single"

    input:
        tuple val(meta_study_id), val(meta_parameters), path(gwas_input)

    output:
        tuple val(meta_study_id), emit: validation

    script:
    """
    required_columns=\"${meta_parameters.rsid_lab} ${meta_parameters.chr_lab} ${meta_parameters.pos_lab} ${meta_parameters.a1_lab} ${meta_parameters.a0_lab} ${meta_parameters.effect_lab} ${meta_parameters.se_lab}\"
    
    \"${projectDir}/bin/validate_columns.sh\" \"\$required_columns\" \"${gwas_input}\"
    """
}
