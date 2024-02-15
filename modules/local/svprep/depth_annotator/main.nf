process SVPREP_DEPTH_ANNOTATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.3--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sv-prep:1.2.3--hdfd78af_1' }"

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
    svprep \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.svprep.depth.DepthAnnotator \\
        ${args} \\
        -input_vcf ${vcf} \\
        -samples ${labels_arg} \\
        -bam_files ${bams_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -threads ${task.cpus} \\
        -output_vcf ${meta.tumor_id}.gridss.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svprep: \$(svprep -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.gridss.vcf.gz
    touch ${meta.tumor_id}.gridss.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
