process AMBER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/amber:4.0.rc1--0'

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    val ref_genome_ver
    path heterozygous_sites
    path target_region_bed

    output:
    tuple val(meta), path('amber/'), emit: amber_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    def target_regions_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -tumor ${meta.tumor_id} \\
            -tumor_bam ${tumor_bam} \\
            ${reference_arg} \\
            ${reference_bam_arg} \\
            ${target_regions_bed_arg} \\
            -ref_genome_version ${ref_genome_ver} \\
            -loci ${heterozygous_sites} \\
            -threads ${task.cpus} \\
            -output_dir amber/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: \$(java -jar ${task.ext.jarPath} -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amber/
    touch amber/placeholder
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
