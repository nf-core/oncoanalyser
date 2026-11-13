process BWAMEM3_ALIGN_RNA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-6bd47f6a3712e01a2e962d032a73db77d068fcbb:7ee0b42c43307d1a57d339331336746a337da55b-0' :
        'biocontainers/mulled-v2-6bd47f6a3712e01a2e962d032a73db77d068fcbb:7ee0b42c43307d1a57d339331336746a337da55b-0' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    path genome_bwamem2_index

    output:
    tuple val(meta), path('*.bam')                                , topic: bwamem3_align_rna_bam
    tuple val(meta), val('bwamem3_align_rna'), path('.command.*') , topic: command_files
    path 'versions.yml'                                           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''

    def output_fn = meta.split ? "${meta.split}.${meta.output_file_id}.bam" : "${meta.output_file_id}.bam"

    // NOTE(LN): TARS consumes read groups:
    // - Output BAM must therefore stay name-grouped, not coordinate sorted
    // - No BAI file created - that would require a coordinate sorted BAM
    """
    ln -fs \$(find -L ${genome_bwamem2_index} -type f) ./

    # NOTE(RJ): compat mode writes the same records as bwa-mem2 and reads the bwa-mem2 index
    bwa-mem3 mem \\
        --compat=bwa-mem2 \\
        ${args} \\
        -Y \\
        -K 100000000 \\
        -R '${meta.rg_line}' \\
        -t ${task.cpus} \\
        -T 19 \\
        -h 75 \\
        \\
        ${genome_fasta} \\
        ${reads_fwd} \\
        ${reads_rev} | \\
        \\
        sambamba view \\
            ${args2} \\
            --sam-input \\
            --format bam \\
            --nthreads ${task.cpus} \\
            --output-filename ${output_fn} \\
            /dev/stdin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem3: \$(bwa-mem3 version 2>&1 | head -n1)
        sambamba:\$(sambamba --version 2>&1 | sed -n '/^sambamba / { s/^.* //p }' | head -n1)
    END_VERSIONS
    """

    stub:
    def output_fn = meta.split ? "${meta.split}.${meta.output_file_id}.bam" : "${meta.output_file_id}.bam"

    """
    touch ${output_fn}

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
