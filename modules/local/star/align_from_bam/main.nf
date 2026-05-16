process STAR_ALIGN_FROM_BAM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-69a6f67cb46e41b4c393f71634a9956d5e31f3e9:46745d95bbdc75d9503849416a66ac6555567ff0-0' :
        'quay.io/biocontainers/mulled-v2-69a6f67cb46e41b4c393f71634a9956d5e31f3e9:46745d95bbdc75d9503849416a66ac6555567ff0-0' }"

    input:
    tuple val(meta), path(bam_input)   // BAM or CRAM to be realigned
    path genome_star_index

    output:
    tuple val(meta), path('Aligned.out.bam'), emit: bam
    path 'versions.yml'          , emit: versions
    path '.command.*'            , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    # Create a collated version of the bam on local disk (annoying use of disk, but it's RNA so likely not so big)
    # Unfortunately there's a FIFO deadlock if this is all piped together
    samtools collate -@ ${task.cpus} -u -o collated.bam ${bam_input}

    STAR \\
        ${args} \\
        --readFilesIn \\
          <(samtools fastq -@ ${task.cpus} -n -F 0x900 -1 /dev/stdout -2 /dev/null -0 /dev/null -s /dev/null collated.bam) \\
          <(samtools fastq -@ ${task.cpus} -n -F 0x900 -1 /dev/null -2 /dev/stdout -0 /dev/null -s /dev/null collated.bam) \\
        --genomeDir ${genome_star_index} \\
        --runThreadN ${task.cpus} \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --alignSplicedMateMapLmin 35 \\
        --alignSplicedMateMapLminOverLmate 0.33 \\
        --chimJunctionOverhangMin 10 \\
        --chimOutType WithinBAM SoftClip \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreMin 1 \\
        --chimScoreSeparation 1 \\
        --chimSegmentMin 10 \\
        --chimSegmentReadGapMax 3 \\
        --limitOutSJcollapsed 3000000 \\
        --outBAMcompression 0 \\
        --outFilterMatchNmin 35 \\
        --outFilterMatchNminOverLread 0.33 \\
        --outFilterMismatchNmax 3 \\
        --outFilterMultimapNmax 10 \\
        --outFilterScoreMinOverLread 0.33 \\
        --outSAMattributes All \\
        --outSAMattrRGline ID:${meta.read_group} SM:${meta.sample_id} \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --runRNGseed 0

    # Wait for background samtools process and propagate any error
    wait \$!

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch Aligned.out.bam

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
