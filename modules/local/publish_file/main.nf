process PUBLISH_FILE {
  tag "${collected_file.fileName.name}"
  label "process_single"

  publishDir path: { destination }, mode: 'copy', overwrite: true

  input:
    path(collected_file, stageAs: 'input/*')
    val(destination)

  output:
    path "${collected_file.fileName.name}", emit: published_file

  script:
    """
    cp -- "${collected_file}" "${collected_file.fileName.name}"
    """
}