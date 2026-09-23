process BWAMEM3_INDEX {
    tag "$fasta"
    label 'process_single'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/675fb6ba74af9bf75f4f7d8819f8f0590d93a625444fc4d8e000b4bddf6d8d4e/data' :
        'community.wave.seqera.io/library/bwa-mem3_samtools_sambamba:630ef504d2934496' }"

    input:
    path fasta
    path alt

    output:
    path 'bwa-mem2_index'                                   , topic: bwamem3_index
    tuple val([:]), val('bwamem3_index'), path('.command.*'), topic: command_files
    path 'versions.yml'                                     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''

    // NOTE: bwa-mem3 reads and writes the bwa-mem2 index format, so the output directory keeps the bwa-mem2_index name
    """
    mkdir -p bwa-mem2_index/
    bwa-mem3 \\
        index \\
        $args \\
        -t ${task.cpus} \\
        -p bwa-mem2_index/${prefix} \\
        $fasta

    # Include ALT file where necessary
    if [[ -n "${alt}" ]]; then
        ln -s ../${alt} bwa-mem2_index/;
    fi;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem3: \$(bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"

    """
    mkdir -p bwa-mem2_index/

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
