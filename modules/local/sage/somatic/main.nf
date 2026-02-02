// NOTE(SW): logic that determines BQR outputs assumes '-output_vcf' is a path that includes at least leading one directory

process SAGE_SOMATIC {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:4.2--hdfd78af_0' :
        'biocontainers/hmftools-sage:4.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(donor_bam), path(tumor_bai), path(normal_bai), path(donor_bai), path(redux_tsvs)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path sage_pon
    path sage_known_hotspots_somatic
    path sage_highconf_regions
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource
    val sequencing_type
    val targeted_mode

    output:
    tuple val(meta), path('somatic/*.sage.somatic.vcf.gz'), path('somatic/*.sage.somatic.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('somatic/')                                                                , emit: sage_dir
    path 'versions.yml'                                                                              , emit: versions
    path '.command.*'                                                                                , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_ids = [meta.normal_id, meta.donor_id].findAll { it }
    def reference_bams = [normal_bam, donor_bam].findAll { it }.collect { it.toString() }

    def reference_arg = reference_ids ? "-reference ${reference_ids.join(',')}" : ''
    def reference_bam_arg = reference_bams ? "-reference_bam ${reference_bams.join(',')}" : ''
    def ref_sample_count_arg = reference_ids ? "-ref_sample_count ${reference_ids.size()}" : ''

    // Tumor in normal contamination (TINC): only for WGS tumor/normal samples
    def run_tinc_arg = ''
    def write_fit_variants_arg = ''
    def gnomad_arg = ''
    def pon_file_arg = ''

    if (!targeted_mode && tumor_bam && normal_bam) {
        run_tinc_arg = '-run_tinc'
        pon_file_arg = "-pon_file ${sage_pon}"
        write_fit_variants_arg = '-write_fit_variants'

        if (genome_ver.toString() == '37') {
            gnomad_arg = "-gnomad_freq_file ${gnomad_resource}"
        } else if (genome_ver.toString() == '38') {
            gnomad_arg = "-gnomad_freq_dir ${gnomad_resource}"
        } else {
            error "got bad genome version: ${genome_ver}"
        }
    }

    // NOTE(SW): use of ternary inexplicitly causes a 'variable already defined in scope error'
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
        -tumor_bam ${tumor_bam} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -hotspots ${sage_known_hotspots_somatic} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -high_confidence_bed ${sage_highconf_regions} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -sequencing_type ${sequencing_type} \\
        -include_mt \\
        ${pon_file_arg} \\
        ${gnomad_arg} \\
        ${run_tinc_arg} \\
        ${high_depth_mode_arg} \\
        ${write_fit_variants_arg} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
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

    ${ (meta.normal_id != null) ? "touch somatic/${meta.normal_id}.sage.bqr.{png,tsv}" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
