process TARS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-tars:1.0--hdfd78af_0' :
        'biocontainers/hmftools-tars:1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(alns)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path contig_sidecar
    path unmap_regions_rna

    output:
    tuple val(meta), path('*.tars.bam'), path('*.tars.bam.bai'), topic: tars_bam
    tuple val(meta), path('*.tars.summary.tsv')                , topic: tars_summary
    tuple val(meta), val('tars'), path('.command.*')           , topic: command_files
    path 'versions.yml'                                        , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    tars \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_bam ${alns.join(',')} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -contig_sidecar ${contig_sidecar} \\
        -rna_unmap_regions ${unmap_regions_rna} \\
        -bamtool \$(which samtools) \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tars: \$(tars -version | sed -n '/^Tars version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.tars.bam
    touch ${meta.sample_id}.tars.bam.bai
    touch ${meta.sample_id}.tars.summary.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
