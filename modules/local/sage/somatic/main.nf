nextflow.enable.types = true

// NOTE(SW): logic that determines BQR outputs assumes '-output_vcf' is a path that includes at least leading one directory

process SAGE_SOMATIC {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:5.0.2--hdfd78af_0' :
        'biocontainers/hmftools-sage:5.0.2--hdfd78af_0' }"

    input:
    tuple(
        meta: Record,
        tumor_aln: Path,
        tumor_bai: Path,
        normal_aln: Path?,
        normal_bai: Path?,
        donor_aln: Path?,
        donor_bai: Path?,
        redux_tsvs: List<Path>,
    )
    genome_fasta: Path
    genome_ver: String
    genome_fai: Path
    genome_dict: Path
    sage_pon: Path
    sage_known_hotspots_somatic: Path
    sage_highconf_regions: Path
    driver_gene_panel: Path
    ensembl_data_resources: Path
    gnomad_resource: Path
    sequencing_platform: String
    targeted_mode: Boolean

    topic:
    tuple(meta, file('somatic/'))                    >> 'sage_somatic_dir'
    tuple(meta, 'sage_somatic', files('.command.*')) >> 'command_files'
    file('versions.yml')                             >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_ids = []
    if (meta.normal_id != null) { reference_ids.add(meta.normal_id) }
    if (meta.donor_id != null) { reference_ids.add(meta.donor_id) }
    def reference_arg = reference_ids.size() > 0 ? "-reference ${reference_ids.join(',')}" : ''
    def ref_sample_count_arg = reference_ids.size() > 0 ? "-ref_sample_count ${reference_ids.size()}" : ''

    def reference_alns = []
    if (normal_aln) { reference_alns.add(normal_aln.toString()) }
    if (donor_aln) { reference_alns.add(donor_aln.toString()) }
    def reference_bam_arg = reference_alns.size() > 0 ? "-reference_bam ${reference_alns.join(',')}" : ''

    def include_mt_arg = targeted_mode ? '' : '-include_mt'

    // Set TINC when if conditions
    def tinc_args = ''

    def should_run_tinc_wgs_tn = ! targeted_mode && tumor_aln && normal_aln
    def should_run_tinc_seq_type = sequencing_platform.toLowerCase() == 'illumina'
    def should_run_tinc = should_run_tinc_wgs_tn && should_run_tinc_seq_type

    if (should_run_tinc) {

        def gnomad_arg
        if (genome_ver == '38') {
            gnomad_arg = "-gnomad_freq_dir ${gnomad_resource}"
        } else {
            gnomad_arg = "-gnomad_freq_file ${gnomad_resource}"
        }

        tinc_args = "-run_tinc -write_fit_variants -pon_file ${sage_pon} ${gnomad_arg}"

    }

    // NOTE(SW): use of ternary inexplicitly causes a 'variable already defined in scope error'
    def high_depth_mode_arg
    if (targeted_mode) {
        high_depth_mode_arg = '-high_depth_mode'
    } else {
        high_depth_mode_arg = ''
    }

    """
    mkdir -p somatic/

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
        -driver_gene_panel ${driver_gene_panel} \\
        -high_confidence_bed ${sage_highconf_regions} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        ${include_mt_arg} \\
        ${tinc_args} \\
        ${high_depth_mode_arg} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_vcf somatic/${meta.tumor_id}.sage.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed -n '/^Sage version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p somatic/

    gzip <<< '' > somatic/${meta.tumor_id}.sage.somatic.vcf.gz
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz.tbi
    touch somatic/${meta.tumor_id}.gene.coverage.tsv
    touch somatic/${meta.tumor_id}.sage.bqr.png
    touch somatic/${meta.tumor_id}.sage.bqr.tsv

    ${ (meta.normal_id != null) ? "touch somatic/${meta.normal_id}.sage.bqr.{png,tsv}" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
