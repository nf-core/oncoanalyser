// NOTE(SW): the --db argument for the virusbreakend command must have a trailing slash if it is a symlink

process VIRUSBREAKEND {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "nf-core/gridss:2.13.2--1"

    input:
    tuple val(meta), path(bam)
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_gridss_index
    path virusbreakenddb
    path gridss_config

    output:
    tuple val(meta), path("*.summary.tsv"), emit: tsv
    path "*.virusbreakend.vcf"            , emit: vcf
    path 'versions.yml'                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    """
    # Symlink indices next to assembly FASTA
    ln -sf \$(find -L ${genome_gridss_index} -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)') ./

    virusbreakend \\
        ${args} \\
        --gridssargs "--jvmheap ${Math.round(task.memory.bytes * xmx_mod)}" \\
        --threads ${task.cpus} \\
        --db ${virusbreakenddb.toString().replaceAll("/\$", "")}/ \\
        --output ${meta.sample_id}.virusbreakend.vcf \\
        --reference ${genome_fasta} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed -n '/-gridss\$/ { s/-gridss//p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.virusbreakend.vcf ${meta.sample_id}.summary.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
