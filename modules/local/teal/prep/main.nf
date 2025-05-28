process TEAL_PREP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-teal:1.3.3--hdfd78af_0' :
        'biocontainers/hmftools-teal:1.3.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    val genome_ver

    output:
    tuple val(meta), path("teal_bam/${meta.tumor_id}.teal.telbam.bam") , emit: tumor_teal_bam
    tuple val(meta), path("teal_bam/${meta.normal_id}.teal.telbam.bam"), emit: normal_teal_bam, optional: true
    path 'versions.yml'                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def tumor_arg = tumor_bam ? "-tumor ${meta.tumor_id}": ''
    def tumor_bam_arg = tumor_bam ? "-tumor_bam ${tumor_bam}": ''

    def reference_arg = normal_bam ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    """
    teal \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.teal.TealPipelineTelbamApp \\
        ${args} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${tumor_arg} \\
        ${tumor_bam_arg} \\
        -ref_genome_version ${genome_ver} \\
        -threads ${task.cpus} \\
        -output_dir teal_bam/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teal: \$(teal -output_dir ./ | sed -n '/ Teal version:/ { s/^.*: //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p teal/

    ${ (meta.tumor_id != null) ? "touch teal/${meta.tumor_id}.teal.telbam.bam" : '' }
    ${ (meta.normal_id != null) ? "touch teal/${meta.normal_id}.teal.telbam.bam" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
