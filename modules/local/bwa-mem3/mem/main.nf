process BWAMEM3_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/675fb6ba74af9bf75f4f7d8819f8f0590d93a625444fc4d8e000b4bddf6d8d4e/data' :
        'community.wave.seqera.io/library/bwa-mem3_samtools_sambamba:630ef504d2934496' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    path genome_bwamem2_index

    output:
    tuple val(meta), path('*.bam'), path('*.bai')            , topic: bwamem3_align_bam
    tuple val(meta), val('bwamem3_align'), path('.command.*'), topic: command_files
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

    bwa-mem3 mem \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem3: \$(bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p')
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
