process MARKDUPS {
    tag "${meta_bam.id}"

    // TODO(MC): Resources required?
    // label 'process_low'

    container 'docker.io/scwatts/markdups:1.1.rc1'

    input:
    tuple val(meta_bam), path(bams), path(bais)
    path genome_fasta
    path genome_fai
    path genome_dict
    path unmap_regions

    output:
    tuple val(meta_bam), path('*bam'), path('*bai'), emit: bam
    path '*.tsv'

    // TODO(MC): Make sure this is in each.
    when:
    task.ext.when == null || task.ext.when

    // TODO(MC): Versions in each.
    // path 'versions.yml'         , emit: versions

    // script:
    // def args = task.ext.args ?: ''

    // // TODO(SW): implement process
    // """
    // echo bar

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     markdups: foo
    // END_VERSIONS
    // """

    // stub:
    // // TODO(SW): implement stub
    // """
    // touch bar
    // echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    // """

    // # TODO(MC): Umi flags
    //     # -multi_bam \\
    //     # -umi_enabled \\
    //     # -umi_duplex \\
    //     # -umi_duplex_delim _ \\
    //     # -umi_base_diff_stats \\

    script:
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
        \\
        -unmap_regions ${unmap_regions} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version 37 \\
        \\
        -write_stats \\
        -threads 16 \\
        \\
        -output_bam ${meta_bam.sample_id}.mark_dups.bam
    """

    stub:
    """
    touch ${meta_bam.sample_id}.mark_dups.bam
    touch ${meta_bam.sample_id}.mark_dups.bam.bai
    touch ${meta_bam.sample_id}.duplicate_freq.tsv
    """

    // # TODO(MC):
    // # touch ${meta_bam.sample_id}.umi_coord_freq.tsv
    // # touch ${meta_bam.sample_id}.umi_edit_distance.tsv
    // # touch ${meta_bam.sample_id}.umi_nucleotide_freq.tsv
}
