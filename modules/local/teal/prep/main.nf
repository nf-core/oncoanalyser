process TEAL_PREP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-teal:1.4--hdfd78af_0' :
        'biocontainers/hmftools-teal:1.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_aln), path(tumor_idx), path(normal_aln), path(normal_idx)
    path genome_fasta
    val genome_ver
    path genome_fai
    val sequencing_platform

    output:
    tuple val(meta), path("teal_bam/${meta.tumor_id}.teal.telbam{.bam,.bam.bai}") , topic: teal_prep_tumor_bam
    tuple val(meta), path("teal_bam/${meta.normal_id}.teal.telbam{.bam,.bam.bai}"), topic: teal_prep_normal_bam, optional: true
    tuple val(meta), val('teal_prep'), path('.command.*')                         , topic: command_files
    path 'versions.yml'                                                           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def tumor_arg = ''
    def tumor_bam_arg = ''
    if (tumor_aln) {
        tumor_arg = "-tumor ${meta.tumor_id}"
        tumor_bam_arg = "-tumor_bam ${tumor_aln}"
    }

    def reference_arg = ''
    def reference_bam_arg = ''
    if (normal_aln) {
        reference_arg = "-reference ${meta.normal_id}"
        reference_bam_arg = "-reference_bam ${normal_aln}"
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
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir teal_bam/

    if [[ -e "${tumor_aln}" ]]; then samtools index teal_bam/${meta.tumor_id}.teal.telbam.bam; fi
    if [[ -e "${normal_aln}" ]]; then samtools index teal_bam/${meta.normal_id}.teal.telbam.bam; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teal: \$(teal -version | sed -n '/Teal version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
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
