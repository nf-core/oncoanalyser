// NOTE(SW): the --db argument for the virusbreakend command must have a trailing slash if it is a symlink

process VIRUSBREAKEND {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "nf-core/gridss:2.13.2--1"

    input:
    tuple val(meta), path(aln)
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_gridss_index
    path virusbreakenddb
    path gridss_config

    output:
    tuple val(meta), path('*.summary.tsv')                   , topic: virusbreakend_tsv
    tuple val(meta), path('*.virusbreakend.vcf')             , topic: virusbreakend_vcf
    tuple val(meta), val('virusbreakend'), path('.command.*'), topic: command_files
    path 'versions.yml'                                      , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    """
    # Symlink indices next to assembly FASTA
    ln -sf \$(find -L ${genome_gridss_index} -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)') ./

    # NOTE(SW): a htslib-compatible reference cache must be provided for CRAMs with stale UR fields
    if [[ "${aln}" == *.cram ]]; then
        seq_cache_populate.pl -root ref_cache/ -subdirs 2 ${genome_fasta}
        export REF_CACHE=ref_cache/%2s/%2s/%s
    fi

    virusbreakend \\
        ${args} \\
        --gridssargs "--jvmheap ${Math.round(task.memory.bytes * xmx_mod)}" \\
        --threads ${task.cpus} \\
        --db ${virusbreakenddb.toString().replaceAll("/\$", "")}/ \\
        --output ${meta.sample_id}.virusbreakend.vcf \\
        --reference ${genome_fasta} \\
        ${aln}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed -n '/-gridss\$/ { s/-gridss//p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
        bcftools: \$(bcftools --version | sed -n '/^bcftools / { s/^.* //p }')
        bwa: \$(bwa 2>&1 | sed -n '/Version/ { s/^Version: //p }')
        repeatmasker: \$(RepeatMasker -v | sed -n '/version/ { s/Repeat.* //p }')
        kraken2: \$(kraken2 -v | sed -n '/version/ { s/^Kraken version //p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        r-structuralvariantannotation: \$(Rscript -e 'packageVersion("StructuralVariantAnnotation") |> as.character() |> writeLines()')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.virusbreakend.vcf ${meta.sample_id}.summary.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
