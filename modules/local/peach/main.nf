process PEACH {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/scwatts/peach:1.6--0'

    input:
    tuple val(meta), path(germline_vcf)
    val genome_ver
    path panel

    output:
    tuple val(meta), path('*.genotype.tsv'), emit: genotype
    path '*.calls.tsv'                     , emit: calls
    path 'versions.yml'                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tool_version = '1.6'

    """
    python ${task.ext.scriptPath} \\
        --sample_t_id ${meta.tumor_id} \\
        --sample_r_id ${meta.normal_id} \\
        --vcf ${germline_vcf} \\
        --vcf_reference_assembly_version V${genome_ver} \\
        --panel ${panel} \\
        --tool_version ${tool_version} \\
        --outputdir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peach: ${tool_version}
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.tumor_id}.peach.calls.tsv"
    touch "${meta.tumor_id}.peach.genotype.tsv"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
