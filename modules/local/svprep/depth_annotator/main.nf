process SVPREP_DEPTH_ANNOTATOR {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/svprep:1.2.1--0'

    input:
    tuple val(meta), path(bams), path(bais), path(vcf), val(labels)
    path genome_fasta
    val genome_ver

    output:
    tuple val(meta), path("${meta.tumor_id}.gridss.vcf.gz"), emit: vcf
    path("${meta.tumor_id}.gridss.vcf.gz.tbi")             , emit: tbi
    path 'versions.yml'                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def labels_list = labels instanceof List ? labels : [labels]
    def labels_arg = labels_list.join(',')
    // NOTE(SW): Nextflow implicitly casts List<TaskPath> to an atomic TaskPath, hence the required check below
    def bams_list = bams instanceof List ? bams : [bams]
    def bams_arg = "${bams_list.join(',')}"

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} com.hartwig.hmftools.svprep.depth.DepthAnnotator \\
            ${args} \\
            -input_vcf ${vcf} \\
            -samples ${labels_arg} \\
            -bam_files ${bams_arg} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -threads ${task.cpus} \\
            -output_vcf ${meta.tumor_id}.gridss.vcf.gz

    # NOTE(SW): partially hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svprep: 1.2.1
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.gridss.vcf.gz
    touch ${meta.tumor_id}.gridss.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
