process CIDER {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cider:1.2--hdfd78af_0' :
        'biocontainers/hmftools-cider:1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(aln), path(idx)
    file genome_fasta
    val genome_ver
    file genome_fai
    file genome_dict
    file genome_img

    output:
    tuple val(meta), path('cider/*')                 , topic: cider_results
    tuple val(meta), val('cider'), path('.command.*'), topic: command_files
    path 'versions.yml'                              , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    cider \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.cider.CiderApplication \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam ${aln} \\
        -ref_genome_version ${genome_ver} \\
        -ref_genome ${genome_fasta} \\
        -write_cider_bam \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir cider/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cider: \$(cider -version | sed -n '/^Cider version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cider/

    touch cider/${meta.sample_id}.cider.bam
    gzip <<< '' > cider/${meta.sample_id}.cider.alignment_match.tsv.gz
    gzip <<< '' > cider/${meta.sample_id}.cider.layout.gz
    touch cider/${meta.sample_id}.cider.locus_stats.tsv
    gzip <<< '' > cider/${meta.sample_id}.cider.vdj.tsv.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
