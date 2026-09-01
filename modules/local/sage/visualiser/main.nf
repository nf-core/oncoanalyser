// NOTE(SW): logic that determines BQR outputs assumes '-output_vcf' is a path that includes at least leading one directory

process SAGE_VISUALISER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:5.0.2--hdfd78af_0' :
        'biocontainers/hmftools-sage:5.0.2--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(tumor_aln),
        path(normal_aln),
        path(donor_aln),
        path(tumor_idx),
        path(normal_idx),
        path(donor_idx),
        path(redux_tsvs),
        path(purple_vcf),
        path(purple_vcf_tbi)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path sage_pon
    path sage_known_hotspots_somatic
    path sage_highconf_regions
    path ensembl_data_resources
    val targeted_mode

    output:
    tuple val(meta), path('sage_vis/')                         , topic: sage_visualiser_dir
    tuple val(meta), val('sage_visualiser'), path('.command.*'), topic: command_files
    path 'versions.yml'                                        , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_ids = []
    if (meta.containsKey('normal_id')) { reference_ids.add(meta.normal_id) }
    if (meta.containsKey('donor_id')) { reference_ids.add(meta.donor_id) }
    def reference_arg = reference_ids.size() > 0 ? "-reference ${reference_ids.join(',')}" : ''
    def ref_sample_count_arg = reference_ids.size() > 0 ? "-ref_sample_count ${reference_ids.size()}" : ''

    def reference_alns = []
    if (normal_aln) { reference_alns.add(normal_aln.toString()) }
    if (donor_aln) { reference_alns.add(donor_aln.toString()) }
    def reference_bam_arg = reference_alns.size() > 0 ? "-reference_bam ${reference_alns.join(',')}" : ''

    def include_mt_arg = targeted_mode ? '' : '-include_mt'

    """
    mkdir -p sage_vis/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${ref_sample_count_arg} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_aln} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -hotspots ${sage_known_hotspots_somatic} \\
        -high_confidence_bed ${sage_highconf_regions} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        ${include_mt_arg} \\
        -vis_purple_vcf ${purple_vcf} \\
        -vis_output_dir sage_vis/ \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_vcf sage_vis/${meta.tumor_id}.sage.vis.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed -n '/^Sage version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sage_vis/

    gzip <<< '' > sage_vis/${meta.tumor_id}.sage.vis.vcf.gz
    touch sage_vis/${meta.tumor_id}.sage.vis.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
