process WRITE_REFERENCE_DATA {
    tag "${fp.name}"

    input:
    path fp
    val workflow_version

    output:
    path fp, includeInputs: true

    script:
    """
    """
}
