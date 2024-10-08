// NOTE(SW): logic that determines BQR outputs assumes '-output_vcf' is a path that includes at least leading one directory

process SAGE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:4.0_beta--hdfd78af_0' :
        'biocontainers/hmftools-sage:4.0_beta--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(donor_bam), path(tumor_bai), path(normal_bai), path(donor_bai)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path sage_known_hotspots_somatic
    path sage_actionable_panel
    path sage_coverage_panel
    path sage_highconf_regions
    path ensembl_data_resources

    output:
    tuple val(meta), path('somatic/*.sage.somatic.vcf.gz'), path('somatic/*.sage.somatic.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('somatic/')                                                                , emit: sage_dir
    path 'versions.yml'                                                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_ids = []
    if (meta.normal_id != null) reference_ids.add(meta.normal_id)
    if (meta.donor_id != null) reference_ids.add(meta.donor_id)
    def reference_arg = reference_ids.size() > 0 ? "-reference ${String.join(",", reference_ids)}" : ""

    def reference_bams = []
    if (normal_bam) reference_bams.add(normal_bam.toString())
    if (donor_bam) reference_bams.add(donor_bam.toString())
    def reference_bam_arg = reference_bams.size() > 0 ? "-reference_bam ${String.join(",", reference_bams)}" : ""

    def ref_sample_count_arg = "-ref_sample_count ${reference_ids.size()}"

    def run_mode = Utils.getEnumFromString(params.mode, Constants.RunMode)
    def high_depth_mode_arg = (run_mode === Constants.RunMode.TARGETED) ? "-high_depth_mode" : ""

    """
    # Get MSI jitter files
    mkdir -p redux/

    symlink_redux_tsvs_from () {
        # TSV files are stored in the same dir as the BAM files
        bam_file=\$(realpath \$1)
        bam_dir=\$(dirname \$bam_file)
        ln -sf \$bam_dir/*.tsv* redux/
    }

    symlink_redux_tsvs_from ${tumor_bam}
    ${ (meta.normal_id != null) ? "symlink_redux_tsvs_from ${normal_bam}" : "" }
    ${ (meta.donor_id != null) ? "symlink_redux_tsvs_from ${donor_bam}" : "" }

    # Run SAGE
    mkdir -p somatic/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${ref_sample_count_arg} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        -jitter_param_dir redux/ \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -hotspots ${sage_known_hotspots_somatic} \\
        -panel_bed ${sage_actionable_panel} \\
        -coverage_bed ${sage_coverage_panel} \\
        -high_confidence_bed ${sage_highconf_regions} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        ${high_depth_mode_arg} \\
        -bqr_write_plot \\
        -threads ${task.cpus} \\
        -output_vcf somatic/${meta.tumor_id}.sage.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p somatic/
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz.tbi
    touch somatic/${meta.tumor_id}.gene.coverage.tsv
    touch somatic/${meta.tumor_id}.sage.bqr.png
    touch somatic/${meta.tumor_id}.sage.bqr.tsv
    touch somatic/${meta.normal_id}.sage.bqr.png
    touch somatic/${meta.normal_id}.sage.bqr.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
