process SVPREP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.3--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sv-prep:1.2.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(junctions)
    path genome_fasta
    val genome_ver
    path sv_blocklist
    path known_fusions
    val write_types

    output:
    tuple val(meta), path("*.sorted.bam")           , emit: bam
    tuple val(meta), path("*.sv_prep.junctions.tsv"), emit: junctions
    path 'versions.yml'                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def write_types_arg = write_types ? "-write_types \'${write_types}\'" : ''
    def existing_juction_file_arg = junctions ? "-existing_junction_file ${junctions}" : ''

    """
    svprep \\
        -Xmx${Math.round(task.memory.bytes * 0.75)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam_file ${bam} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -blacklist_bed ${sv_blocklist} \\
        -known_fusion_bed ${known_fusions} \\
        ${write_types_arg} \\
        ${existing_juction_file_arg} \\
        -threads ${task.cpus} \\
        -output_dir ./

    samtools sort \\
        -@ ${task.cpus} \\
        -T ${meta.sample_id}.sv_prep.tmp \\
        -o ${meta.sample_id}.sv_prep.sorted.bam \\
        ${meta.sample_id}.sv_prep.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svprep: \$(svprep -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.sample_id}.sv_prep.sorted.bam"
    touch "${meta.sample_id}.sv_prep.junctions.tsv"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
