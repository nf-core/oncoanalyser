process BWAMEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:544519c4a0ff7e9616a3b44afde1f143c52f10c3-0' :
        'quay.io/biocontainers/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:544519c4a0ff7e9616a3b44afde1f143c52f10c3-0' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    path genome_bwa_index
    path genome_bwa_index_bseq
    path genome_bwa_index_biidx

    output:
    tuple val(meta), path('*.bam'), path('*.bai'), emit: bam
    path 'versions.yml'                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def read_group_tag = "@RG\\tID:${meta.read_group}\\tSM:${meta.sample_id}"
    def output_fn = meta.split ? "${meta.split}.${meta.sample_id}.${meta.read_group}.bam" : "${meta.sample_id}.${meta.read_group}.bam"

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
            --out ${output_fn} \\
            /dev/stdin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null)
        sambamba: \$(sambamba --version 2>&1 | egrep '^sambamba' | head -n 1 | awk '{ print \$NF }')
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
