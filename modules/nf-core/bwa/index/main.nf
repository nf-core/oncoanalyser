process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'biocontainers/bwa:0.7.17--hed695b0_7' }"

    input:
    path fasta
    path alt

    output:
    path bwa_index     , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta.name}"
    def args   = task.ext.args ?: ''

    """
    mkdir -p bwa_index/
    bwa \\
        index \\
        $args \\
        -p bwa_index/${prefix} \\
        $fasta

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -sf ../${alt} bwa_index/;
    fi;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.name}"
    """
    mkdir -p bwa_index/

    touch bwa_index/${prefix}.amb
    touch bwa_index/${prefix}.ann
    touch bwa_index/${prefix}.bwt
    touch bwa_index/${prefix}.pac
    touch bwa_index/${prefix}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
