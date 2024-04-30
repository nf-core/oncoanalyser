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
    path bwamem2       , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''

    """
    mkdir -p bwamem2/
    bwa-mem2 \\
        index \\
        $args \\
        $fasta -p bwamem2/${prefix}

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwamem2/;
    fi;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"

    """
    mkdir -p bwamem2/
    touch bwamem2/${prefix}.0123
    touch bwamem2/${prefix}.ann
    touch bwamem2/${prefix}.pac
    touch bwamem2/${prefix}.amb
    touch bwamem2/${prefix}.bwt.2bit.64

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwamem2/;
    fi;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
