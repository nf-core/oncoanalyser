process CIDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cider:1.0.4--hdfd78af_0' :
        'biocontainers/hmftools-cider:1.0.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val genome_ver
    file human_blastdb

    output:
    tuple val(meta), path('cider/*'), emit: cider_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    cider \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.cider.CiderApplication \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam ${bam} \\
        -blast \$(which blastn | sed 's#/bin/blastn##') \\
        -blast_db ${human_blastdb} \\
        -ref_genome_version ${genome_ver} \\
        -threads ${task.cpus} \\
        -write_cider_bam \\
        -output_dir cider/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cider: \$(cider -ref_genome_version 38 -output_dir ./ | sed -n '/ Cider version: / { s/^.*version: \\([0-9.]\\+\\),.*\$/\\1/p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cider/

    touch cider/${meta.sample_id}.cider.bam
    touch cider/${meta.sample_id}.cider.blastn_match.tsv.gz
    touch cider/${meta.sample_id}.cider.layout.gz
    touch cider/${meta.sample_id}.cider.locus_stats.tsv
    touch cider/${meta.sample_id}.cider.vdj.tsv.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
