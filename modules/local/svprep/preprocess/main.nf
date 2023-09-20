process GRIDSS_PREPROCESS {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/svprep:1.2.2--0'

    input:
    tuple val(meta), path(bam), path(bam_filtered)
    path gridss_config
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_bwa_index_dir, stageAs: 'bwa_index'
    path genome_bwa_index_image
    path genome_gridss_index

    output:
    tuple val(meta), path("gridss_preprocess/${meta.id}.sv_prep.sorted.bam.gridss.working/"), emit: preprocess_dir
    path 'versions.yml'                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config_arg = gridss_config ? "--configuration ${gridss_config}" : ''

    """
    # Symlink BWA indices next to assembly FASTA
    ln -s \$(find -L ${genome_bwa_index_dir} -type f) ./

    gridss_svprep \\
        ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --jar ${task.ext.jarPath} \\
        --steps preprocess \\
        --reference ${genome_fasta} \\
        --workingdir gridss_preprocess/ \\
        --threads ${task.cpus} \\
        ${config_arg} \\
        --labels ${meta.id} \\
        --bams ${bam} \\
        --filtered_bams ${bam_filtered} \\
        --output placeholder

    # NOTE(SW): partially hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(java -cp ${task.ext.jarPath} gridss.CallVariants --version 2>&1 | sed 's/-gridss//')
        svprep: 1.2.1
    END_VERSIONS
    """

    stub:
    """
    mkdir -p gridss_preprocess/${meta.id}.sv_prep.sorted.bam.gridss.working/
    touch gridss_preprocess/${meta.id}.sv_prep.sorted.bam.gridss.working/placeholder
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
