process TEAL_PREP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-teal:1.3.5--hdfd78af_0' :
        'biocontainers/hmftools-teal:1.3.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    val genome_ver

    output:
    tuple val(meta), path("teal_bam/${meta.tumor_id}.teal.telbam{.bam,.bam.bai}") , emit: tumor_bam
    tuple val(meta), path("teal_bam/${meta.normal_id}.teal.telbam{.bam,.bam.bai}"), emit: normal_bam, optional: true
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def tumor_arg = ""
    def tumor_bam_arg = ""
    def tumor_bam_index_command = ""

    if(tumor_bam) {
        tumor_arg = "-tumor ${meta.tumor_id}"
        tumor_bam_arg = "-tumor_bam ${tumor_bam}"
        tumor_bam_index_command = "samtools index teal_bam/${meta.tumor_id}.teal.telbam.bam"
    }

    def reference_arg = ""
    def reference_bam_arg = ""
    def reference_bam_index_command = ""

    if(normal_bam) {
        reference_arg = "-reference ${meta.normal_id}"
        reference_bam_arg = "-reference_bam ${normal_bam}"
        reference_bam_index_command = "samtools index teal_bam/${meta.normal_id}.teal.telbam.bam"
    }

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

    ${tumor_bam_index_command}
    ${reference_bam_index_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teal: \$(teal -version | sed -n '/Teal version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p teal_bam/

    ${ (meta.tumor_id != null) ? "touch teal_bam/${meta.tumor_id}.teal.telbam{.bam,.bam.bai}" : '' }
    ${ (meta.normal_id != null) ? "touch teal_bam/${meta.normal_id}.teal.telbam{.bam,.bam.bai}" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
