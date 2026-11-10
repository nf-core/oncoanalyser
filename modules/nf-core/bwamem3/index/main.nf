process BWAMEM3_INDEX {
    tag "${meta.id}"
    // NOTE bwa-mem3 builds an FM-index with libsais; peak memory scales with the reference size.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e6/e6ad16b0d97b68ea5ec5c9a29c372ffe73f54488bb262a03f8496868488463ad/data'
        : 'community.wave.seqera.io/library/bwa-mem3:0.12.0--d4a9f23b4fb7bcd6'}"

    input:
    tuple val(meta), path(fasta)
    path alt

    output:
    tuple val(meta), path("bwamem3"), emit: index
    tuple val(meta), val('bwamem3_index'), path('.command.*'), topic: command_files
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir bwamem3
    bwa-mem3 \\
        index \\
        $args \\
        -t ${task.cpus} \\
        -p bwamem3/${prefix} \\
        $fasta

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwamem3/;
    fi;
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwamem3
    touch bwamem3/${prefix}.0123
    touch bwamem3/${prefix}.amb
    touch bwamem3/${prefix}.ann
    touch bwamem3/${prefix}.bwt.2bit.64
    touch bwamem3/${prefix}.pac

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwamem3/;
    fi;
    """
}
