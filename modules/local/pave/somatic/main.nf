process PAVE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/pave:1.4.3--0'

    input:
    tuple val(meta), path(sage_vcf)
    path genome_fasta
    path genome_fai
    val genome_ver
    path sage_pon
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pon_filters
    def gnomad_args
    if (genome_ver == '37') {
        pon_filters = 'HOTSPOT:10:5;PANEL:6:5;UNKNOWN:6:0'
        gnomad_args = "-gnomad_freq_file ${gnomad_resource}"
    } else if (genome_ver == '38') {
        pon_filters = 'HOTSPOT:5:5;PANEL:2:5;UNKNOWN:2:0'
        gnomad_args = "-gnomad_freq_dir ${gnomad_resource} -gnomad_load_chr_on_demand"
    } else {
        log.error "got bad genome version: ${genome_ver}"
        System.exit(1)
    }

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            -sample ${meta.id} \\
            -vcf_file ${sage_vcf} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -pon_file ${sage_pon} \\
            -pon_filters "${pon_filters}" \\
            -driver_gene_panel ${driver_gene_panel} \\
            -mappability_bed ${segment_mappability} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            ${gnomad_args} \\
            -read_pass_only \\
            -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(java -jar ${task.ext.jarPath} 2>&1 | sed -n 's/^.*version: //p')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sage.pave_somatic.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
