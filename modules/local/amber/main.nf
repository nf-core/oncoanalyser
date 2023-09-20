process AMBER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/amber:3.9--3'

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    val ref_genome_ver
    path heterozygous_sites

    output:
    tuple val(meta), path('amber/'), emit: amber_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -tumor ${meta.tumor_id} \\
            -tumor_bam ${tumor_bam} \\
            ${reference_arg} \\
            ${reference_bam_arg} \\
            -ref_genome_version ${ref_genome_ver} \\
            -loci ${heterozygous_sites} \\
            -threads ${task.cpus} \\
            -output_dir amber/

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: 3.9
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amber/
    touch amber/placeholder
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
