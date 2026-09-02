process BWAMEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:596c0d6a494faa218562f2be03af2714d454da4f-0' :
        'biocontainers/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:596c0d6a494faa218562f2be03af2714d454da4f-0' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    path genome_bwamem2_index

    output:
    tuple val(meta), path('*.bam'), path('*.bai')            , topic: bwamem2_align_bam
    tuple val(meta), val('bwamem2_align'), path('.command.*'), topic: command_files
    path 'versions.yml'                                      , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''

    def output_fn = meta.split ? "${meta.split}.${meta.output_file_id}.bam" : "${meta.output_file_id}.bam"

    """
    ln -fs \$(find -L ${genome_bwamem2_index} -type f) ./

    bwa-mem2 mem \\
        ${args} \\
        -Y \\
        -K 100000000 \\
        -R '${meta.rg_line}' \\
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

    # NOTE(SW): bwa-mem2 version hardcoded as 2.3 reports the wrong version, see https://github.com/bwa-mem2/bwa-mem2/issues/276
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: 2.3
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
        sambamba: \$(sambamba --version 2>&1 | sed -n '/^sambamba / { s/^.* //p }' | head -n1)
    END_VERSIONS
    """

    stub:
    def output_fn = meta.split ? "${meta.split}.${meta.output_file_id}.bam" : "${meta.output_file_id}.bam"

    """
    touch ${output_fn}
    touch ${output_fn}.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
