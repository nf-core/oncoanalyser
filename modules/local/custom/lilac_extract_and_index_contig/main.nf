process CUSTOM_EXTRACTCONTIG {
    tag "${contig_name}"
    label 'process_single'

    conda "bwa=0.7.17 samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4908f45b5b676e6a2fed6b1977d445b16b7b8dee:2411d64e0a784487e81828123aaf68a549531e5c-0' :
        'quay.io/biocontainers/mulled-v2-4908f45b5b676e6a2fed6b1977d445b16b7b8dee:2411d64e0a784487e81828123aaf68a549531e5c-0' }"

    input:
    val contig_name
    path genome_fasta
    path genome_fai
    val run

    output:
    path "*extracted.fa"  , emit: contig
    path "*extracted.fa.*", emit: bwa_indices
    path 'versions.yml'   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    samtools faidx \\
        -o ${contig_name}_extracted.fa \\
        ${genome_fasta} \\
        ${contig_name}
    bwa index ${contig_name}_extracted.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${contig_name}_extracted.fa ${contig_name}_extracted.fa.amb
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
