process BAMTOOLS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-bam-tools:1.3--hdfd78af_0' :
        'biocontainers/hmftools-bam-tools:1.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path genome_fasta
    val genome_ver
    path driver_gene_panel
    path ensembl_data_resources

    output:
    tuple val(meta), path("${meta.id}_bamtools/"), emit: metrics_dir
    path 'versions.yml'                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def output_dir = "${meta.id}_bamtools/"

    def log_file_id = "${task.process.split(':')[-1]}.${meta.sample_id}"

    """
    mkdir -p ${output_dir}

    bamtools \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.bamtools.metrics.BamMetrics \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam_file ${bam} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -threads ${task.cpus} \\
        -log_level ${params.module_log_level} \\
        -output_dir ${output_dir}

    for log_file_ext in sh out err; do
        cp .command.\${log_file_ext} ${output_dir}/command.${log_file_id}.\${log_file_ext}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$(bamtools -version | sed -n '/^BamTools version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_bamtools/

    touch ${meta.id}_bamtools/${meta.sample_id}.bam_metric.summary.tsv;
    touch ${meta.id}_bamtools/${meta.sample_id}.bam_metric.coverage.tsv;
    touch ${meta.id}_bamtools/${meta.sample_id}.bam_metric.frag_length.tsv;
    touch ${meta.id}_bamtools/${meta.sample_id}.bam_metric.flag_counts.tsv;
    touch ${meta.id}_bamtools/${meta.sample_id}.bam_metric.partition_stats.tsv;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
