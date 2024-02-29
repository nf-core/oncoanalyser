process BWAMEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    // TODO(SW): create BioContainers multi-package image when appropriate
    container 'docker.io/scwatts/bwa-mem2:2.2.1'

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    path genome_bwa_index
    path genome_bwa_index_bseq
    path genome_bwa_index_biidx

    output:
    tuple val(meta), path('*.bam'), emit: bam
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def read_group_tag = "@RG\\tID:${meta.read_group}\\tSM:${meta.sample_id}"

    """
    ln -fs \$(find -L ${genome_bwa_index} -type f) ./

    bwa-mem2 mem \\
        -Y \\
        -R '${read_group_tag}' \\
        -t ${task.cpus} \\
        ${genome_fasta} \\
        ${reads_fwd} \\
        ${reads_rev} | \\
        \\
        sambamba view \\
            --sam-input \\
            --format bam \\
            --compression-level 0 \\
            --nthreads ${task.cpus} \\
            /dev/stdin | \\
        \\
        sambamba sort \\
            --nthreads ${task.cpus} \\
            --out ${meta.split}.${meta.sample_id}.${meta.read_group}.bam \\
            /dev/stdin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null)
        sambamba: \$(sambamba --version 2>&1 | egrep '^sambamba' | head -n 1 | awk '{ print \$NF }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.split}.${meta.sample_id}.${meta.read_group}.bam

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
