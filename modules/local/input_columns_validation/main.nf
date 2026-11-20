process INPUT_COLUMNS_VALIDATION {
    label "process_single"

    input:
        tuple val(records), path(gwas_files, stageAs: 'gwas???/*'), path(sdY_files, stageAs: 'sdy???/*')
 
    output:
        val(records), emit: validated_rows
    
    script:
    // Converts each item in records into a JSON string, then joins those strings with newline characters, creating a JSON Lines manifest.
    def manifest_jsonl = records.collect { groovy.json.JsonOutput.toJson(it) }.join('\n')
    // Encodes the manifest text as UTF-8 bytes, then Base64-encodes those bytes and converts the result back to a string.
    def manifest_b64 = manifest_jsonl.getBytes('UTF-8').encodeBase64().toString()
    def args = task.ext.args ?: ''
    """
    base64 --decode > validation_manifest.jsonl <<'MANIFEST_B64'
    ${manifest_b64}
    MANIFEST_B64
    validate_columns.py --manifest validation_manifest.jsonl ${args}
    """

    stub:
    """
    echo "Validation stub for ${records.size()} input rows"
    """
}
