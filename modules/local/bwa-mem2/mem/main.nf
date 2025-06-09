process BWAMEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-bwa-plus:1.0.0--h077b44d_0' :
        'biocontainers/hmftools-bwa-plus:1.0.0--h077b44d_0' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    path genome_bwamem2_index

    output:
    tuple val(meta), path('*.bam'), path('*.bai'), emit: bam
    path 'versions.yml'                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''

    def read_group_tag = "@RG\\tID:${meta.read_group}\\tSM:${meta.sample_id}"
    def output_fn = meta.split ? "${meta.split}.${meta.sample_id}.${meta.read_group}.bam" : "${meta.sample_id}.${meta.read_group}.bam"

    """
    ln -fs \$(find -L ${genome_bwamem2_index} -type f) ./

    bwa-plus mem \\
        ${args} \\
        -Y \\
        -K 100000000 \\
        -R '${read_group_tag}' \\
        -t ${task.cpus} \\
        ${genome_fasta} \\
        ${reads_fwd} \\
        ${reads_rev} | \\
        \\
        sambamba view \\
            ${args2} \\
            --sam-input \\
            --format bam \\
            --compression-level 0 \\
            --nthreads ${task.cpus} \\
            /dev/stdin | \\
        \\
        sambamba sort \\
            ${args3} \\
            --nthreads ${task.cpus} \\
            --out ${output_fn} \\
            /dev/stdin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-plus: \$(bwa-plus version 2>/dev/null)
        sambamba: \$(sambamba --version 2>&1 | sed -n '/^sambamba / { s/^.* //p }' | head -n1)
    END_VERSIONS
    """

    stub:
    def output_fn = meta.split ? "${meta.split}.${meta.sample_id}.${meta.read_group}.bam" : "${meta.sample_id}.${meta.read_group}.bam"

    """
    touch ${output_fn}
    touch ${output_fn}.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
