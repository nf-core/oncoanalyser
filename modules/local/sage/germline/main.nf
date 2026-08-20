nextflow.enable.types = true

process SAGE_GERMLINE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:5.0.2--hdfd78af_0' :
        'biocontainers/hmftools-sage:5.0.2--hdfd78af_0' }"

    input:
    tuple(meta: Map, tumor_aln: Path, tumor_idx: Path, normal_aln: Path?, normal_idx: Path?, redux_tsvs: List<Path>)
    genome_fasta: Path
    genome_ver: String
    genome_fai: Path
    genome_dict: Path
    sage_known_hotspots_germline: Path
    sage_highconf_regions: Path
    driver_gene_panel: Path
    ensembl_data_resources: Path
    sequencing_platform: String
    targeted_mode: Boolean

    topic:
    tuple(meta, file('germline/')) >> 'sage_germline_dir'
    tuple(meta, 'sage_germline', files('.command.*')) >> 'command_files'
    file('versions.yml') >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def high_depth_mode_arg = targeted_mode ? '-high_depth_mode' : ''

    """
    mkdir -p germline/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.normal_id} \\
        -tumor_bam ${normal_aln} \\
        -reference ${meta.tumor_id} \\
        -reference_bam ${tumor_aln} \\
        -ref_sample_count 0 \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -hotspots ${sage_known_hotspots_germline} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -high_confidence_bed ${sage_highconf_regions} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -germline \\
        -panel_only \\
        ${high_depth_mode_arg} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_vcf germline/${meta.tumor_id}.sage.germline.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed -n '/^Sage version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p germline/

    gzip <<< '' > germline/${meta.tumor_id}.sage.germline.vcf.gz
    touch germline/${meta.tumor_id}.sage.germline.vcf.gz.tbi
    touch germline/${meta.tumor_id}.sage.bqr.png
    touch germline/${meta.tumor_id}.sage.bqr.tsv
    touch germline/${meta.normal_id}.sage.bqr.png
    touch germline/${meta.normal_id}.sage.bqr.tsv
    touch germline/${meta.normal_id}.gene.coverage.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
