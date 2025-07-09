process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.19--h577a1d6_0' :
        'biocontainers/bwa:0.7.19--h577a1d6_0' }"

    input:
    path fasta
    path alt

    output:
    path bwa_index     , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''

    def prefix = task.ext.prefix ?: "${fasta.name}"

    """
    mkdir -p bwa_index/
    bwa \\
        index \\
        $args \\
        -p bwa_index/${prefix} \\
        $fasta

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwa_index/;
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

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
