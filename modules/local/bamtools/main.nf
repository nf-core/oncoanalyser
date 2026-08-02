process BAMTOOLS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-bam-tools:1.6.1--hdfd78af_0' :
        'biocontainers/hmftools-bam-tools:1.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(aln), path(idx)
    path genome_fasta
    val genome_ver
    path driver_gene_panel
    path ensembl_data_resources
    path target_regions_bed

    output:
    tuple val(meta), path("bamtools_${meta.sample_id}/"), emit: bamtools_dir
    path 'versions.yml'                                 , emit: versions
    path '.command.*'                                   , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def target_regions_bed_arg = target_regions_bed ? "-regions_file ${target_regions_bed}" : ''

    """
    mkdir -p bamtools_${meta.sample_id}/

    bamtools \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.bamtools.metrics.BamMetrics \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam_file ${aln} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        ${target_regions_bed_arg} \\
        ${log_level_arg} \\
        -threads ${task.cpus} \\
        -output_dir bamtools_${meta.sample_id}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$(bamtools -version | sed -n '/^BamTools version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bamtools_${meta.sample_id}/

    touch bamtools_${meta.sample_id}/${meta.sample_id}.bam_metric.summary.tsv;
    touch bamtools_${meta.sample_id}/${meta.sample_id}.bam_metric.coverage.tsv;
    touch bamtools_${meta.sample_id}/${meta.sample_id}.bam_metric.frag_length.tsv;
    touch bamtools_${meta.sample_id}/${meta.sample_id}.bam_metric.flag_counts.tsv;
    touch bamtools_${meta.sample_id}/${meta.sample_id}.bam_metric.partition_stats.tsv;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
