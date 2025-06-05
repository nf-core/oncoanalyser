process AMBER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-amber:4.1.1--hdfd78af_0' :
        'biocontainers/hmftools-amber:4.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(donor_bam), path(tumor_bai), path(normal_bai), path(donor_bai)
    val genome_ver
    path heterozygous_sites
    path target_region_bed

    output:
    tuple val(meta), path('amber/'), emit: amber_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def reference_ids = []
    if (meta.normal_id != null) reference_ids.add(meta.normal_id)
    if (meta.donor_id != null) reference_ids.add(meta.donor_id)
    def reference_arg = reference_ids.size() > 0 ? "-reference ${String.join(",", reference_ids)}" : ''

    def reference_bams = []
    if (normal_bam) reference_bams.add(normal_bam.toString())
    if (donor_bam) reference_bams.add(donor_bam.toString())
    def reference_bam_arg = reference_bams.size() > 0 ? "-reference_bam ${String.join(",", reference_bams)}" : ''

    def target_regions_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''

    """
    amber \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${target_regions_bed_arg} \\
        -ref_genome_version ${genome_ver} \\
        -loci ${heterozygous_sites} \\
        -threads ${task.cpus} \\
        -output_dir amber/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: \$(amber -version | sed -n '/^Amber version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amber/
    touch amber/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
