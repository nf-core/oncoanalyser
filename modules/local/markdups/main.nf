process MARKDUPS {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-mark-dups:1.1.7--hdfd78af_0' :
        'biocontainers/hmftools-mark-dups:1.1.7--hdfd78af_0' }"

    input:
    tuple val(meta), path(bams), path(bais)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path unmap_regions
    val has_umis

    output:
    tuple val(meta), path('*bam'), path('*bai'), emit: bam
    path 'versions.yml'                        , emit: versions
    path '*.tsv'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def umi_flags = has_umis ? '-umi_enabled -umi_duplex -umi_duplex_delim +' : ''

    """
    markdups \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        \\
        -samtools \$(which samtools) \\
        -sambamba \$(which sambamba) \\
        \\
        -sample ${meta.sample_id} \\
        -input_bam ${bams.join(',')} \\
        \\
        -form_consensus \\
        ${umi_flags} \\
        \\
        -unmap_regions ${unmap_regions} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        \\
        -write_stats \\
        -threads ${task.cpus} \\
        \\
        -output_bam ${meta.sample_id}.markdups.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        markdups: \$(markdups -version | awk '{ print \$NF }')
        sambamba: \$(sambamba --version 2>&1 | egrep '^sambamba' | head -n 1 | awk '{ print \$NF }')
        samtools: \$(samtools --version 2>&1 | egrep '^samtools\\s' | head -n 1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.markdups.bam
    touch ${meta.sample_id}.markdups.bam.bai
    touch ${meta.sample_id}.duplicate_freq.tsv

    if [[ -n "${has_umis}" ]]; then
        touch ${meta.sample_id}.umi_coord_freq.tsv
        touch ${meta.sample_id}.umi_edit_distance.tsv
        touch ${meta.sample_id}.umi_nucleotide_freq.tsv
    fi;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
