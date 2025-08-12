process BWAMEM2_INDEX {
    tag "$fasta"
    label 'process_single'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--he513fc3_0' :
        'biocontainers/bwa-mem2:2.2.1--he513fc3_0' }"

    input:
    path fasta
    path alt

    output:
    path "bwa-mem2_index", emit: index
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''

    """
    mkdir -p bwa-mem2_index/
    bwa-mem2 \\
        index \\
        $args \\
        $fasta -p bwa-mem2_index/${prefix}

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwa-mem2_index/;
    fi;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"

    """
    mkdir -p bwa-mem2_index/
    touch bwa-mem2_index/${prefix}.0123
    touch bwa-mem2_index/${prefix}.ann
    touch bwa-mem2_index/${prefix}.pac
    touch bwa-mem2_index/${prefix}.amb
    touch bwa-mem2_index/${prefix}.bwt.2bit.64

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwa-mem2_index/;
    fi;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
