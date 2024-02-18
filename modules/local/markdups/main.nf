process MARKDUPS {
    tag "${meta_bam.subject_id}__${meta_bam.sample_id}"

    container 'docker.io/scwatts/markdups:1.1.rc1'

    input:
    tuple val(meta_bam), path(bams), path(bais)
    path genome_fasta
    path genome_fai
    path genome_dict
    path unmap_regions
    val has_umis

    output:
    tuple val(meta_bam), path('*bam'), path('*bai'), emit: bam
    path 'versions.yml'                            , emit: versions
    path '*.tsv'

    when:
    task.ext.when == null || task.ext.when

    script:
    def umi_flags = has_umis ? '-umi_enabled -umi_duplex -umi_duplex_delim _ -umi_base_diff_stats' : ''

    // TODO(MC): Ref genome version.

    """
    java \\
      -Xmx${Math.round(task.memory.bytes * 0.95)} \\
      -jar /opt/markdups/markdups.jar \\
        \\
        -samtools \$(which samtools) \\
        -sambamba \$(which sambamba) \\
        \\
        -sample ${meta_bam.sample_id} \\
        -input_bam ${bams.join(',')} \\
        \\
        -form_consensus \\
        -multi_bam \\
        ${umi_flags} \\
        \\
        -unmap_regions ${unmap_regions} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version 37 \\
        \\
        -write_stats \\
        -threads 16 \\
        \\
        -output_bam ${meta_bam.sample_id}.markdups.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: 1.0
        samtools: 1.17
        openjdk: 8
        mark-dups: 1.1
    END_VERSIONS
    """

    stub:
    def umi_output_files = has_umis ? 'touch ${meta_bam.sample_id}.umi_coord_freq.tsv;' +
        ' touch ${meta_bam.sample_id}.umi_edit_distance.tsv;' +
        ' touch ${meta_bam.sample_id}.umi_nucleotide_freq.tsv' : ''

    """
    touch ${meta_bam.sample_id}.markdups.bam
    touch ${meta_bam.sample_id}.markdups.bam.bai
    touch ${meta_bam.sample_id}.duplicate_freq.tsv
    ${umi_output_files}

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
