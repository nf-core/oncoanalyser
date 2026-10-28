process AMBER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-amber:4.3--hdfd78af_0' :
        'biocontainers/hmftools-amber:4.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_aln), path(tumor_idx), path(normal_aln), path(normal_idx), path(donor_aln), path(donor_idx)
    path genome_fasta
    val genome_ver
    path genome_fai
    path heterozygous_sites
    path target_regions_bed
    val tumor_min_depth
    val sequencing_platform

    output:
    tuple val(meta), path('amber/')                  , topic: amber_dir
    tuple val(meta), val('amber'), path('.command.*'), topic: command_files
    path 'versions.yml'                              , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_ids = []
    if (meta.containsKey('normal_id')) { reference_ids.add(meta.normal_id) }
    if (meta.containsKey('donor_id')) { reference_ids.add(meta.donor_id) }
    def reference_arg = reference_ids.size() > 0 ? "-reference ${reference_ids.join(',')}" : ''

    def reference_alns = []
    if (normal_aln) { reference_alns.add(normal_aln.toString()) }
    if (donor_aln) { reference_alns.add(donor_aln.toString()) }

    def reference_bam_arg = reference_alns.size() > 0 ? "-reference_bam ${reference_alns.join(',')}" : ''

    def target_regions_bed_arg = target_regions_bed ? "-target_regions_bed ${target_regions_bed}" : ''

    def tumor_min_depth_arg = tumor_min_depth ? "-tumor_min_depth ${tumor_min_depth}" : ''

    """
    amber \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_aln} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        ${target_regions_bed_arg} \\
        -loci ${heterozygous_sites} \\
        ${tumor_min_depth_arg} \\
        ${log_level_arg} \\
        -threads ${task.cpus} \\
        -output_dir amber/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: \$(amber -version | sed -n '/^Amber version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        bioconductor-copynumber: \$(Rscript -e 'packageVersion("copynumber") |> as.character() |> writeLines()')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amber/

    touch amber/.stub

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
