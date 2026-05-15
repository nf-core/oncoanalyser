process BWAMEM2_ALIGN_FROM_BAM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:596c0d6a494faa218562f2be03af2714d454da4f-0' :
        'biocontainers/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:596c0d6a494faa218562f2be03af2714d454da4f-0' }"

    input:
    tuple val(meta), path(bam_input)        // BAM or CRAM to be realigned
    path genome_fasta
    path genome_bwamem2_index

    output:
    tuple val(meta), path('*.bam'), path('*.bai'), emit: bam
    path 'versions.yml'                           , emit: versions
    path '.command.*'                             , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''

    def read_group_tag = "@RG\\tID:${meta.read_group}\\tSM:${meta.sample_id}"
    def output_fn = "${meta.sample_id}.${meta.read_group}.bam"

    """
    ln -fs \$(find -L ${genome_bwamem2_index} -type f) ./

    samtools collate \\
        -O \\
        -@ ${task.cpus} \\
        ${bam_input} | \\
    samtools fastq \\
        -@ ${task.cpus} \\
        -n \\
        -F 0x900 \\
        -1 /dev/stdout \\
        -2 /dev/stdout \\
        -0 /dev/null \\
        -s /dev/null \\
        - | \\
    \\
    bwa-mem2 mem \\
        ${args} \\
        -Y \\
        -K 100000000 \\
        -p \\
        -R '${read_group_tag}' \\
        -t ${task.cpus} \\
        ${genome_fasta} \\
        /dev/stdin | \\
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

    # Force non-empty output check (helps catch upstream empty input quickly)
    test -s ${output_fn}
    samtools quickcheck ${output_fn}
    samtools index -@ ${task.cpus} ${output_fn}

    # NOTE(SW): bwa-mem2 version hardcoded as 2.3 reports the wrong version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: 2.3
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        sambamba: \$(sambamba --version 2>&1 | sed -n '/^sambamba / { s/^.* //p }' | head -n1)
    END_VERSIONS
    """

    stub:
    def output_fn = "${meta.sample_id}.${meta.read_group}.bam"
    """
    touch ${output_fn}
    touch ${output_fn}.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
