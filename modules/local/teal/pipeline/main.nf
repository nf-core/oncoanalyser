nextflow.enable.types = true

include { TumorNormalMeta } from '../../../../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/records'

process TEAL_PIPELINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-teal:1.4--hdfd78af_0' :
        'biocontainers/hmftools-teal:1.4--hdfd78af_0' }"

    input:
    tuple(
        meta: TumorNormalMeta,
        teal_bam_tumor: Path?,
        tumor_bai_tumor: Path?,
        teal_bam_normal: Path?,
        teal_bai_normal: Path?,
        bamtools_dir_tumor: Path?,
        bamtools_dir_normal: Path?,
        cobalt_dir: Path,
        purple_dir: Path?,
    )
    genome_ver: String
    sequencing_platform: String

    topic:
    tuple(meta, files('teal/*.tsv*'))                 >> 'teal_tsvs'
    tuple(meta, 'teal_pipeline', files('.command.*')) >> 'command_files'
    file('versions.yml')                              >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def tumor_arg = teal_bam_tumor ? "-tumor ${meta.tumor_id}": ''
    def tumor_bam_arg = teal_bam_tumor ? "-tumor_bam ${teal_bam_tumor}": ''
    def tumor_wgs_metrics_arg = bamtools_dir_tumor ? "-tumor_wgs_metrics ${bamtools_dir_tumor}/${meta.tumor_id}.bam_metric.summary.tsv": ''
    def purple_arg = purple_dir ? "-purple_dir ${purple_dir}" : ''

    def reference_arg = teal_bam_normal ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = teal_bam_normal ? "-reference_bam ${teal_bam_normal}" : ''
    def reference_wgs_metrics_arg = bamtools_dir_normal ? "-reference_wgs_metrics ${bamtools_dir_normal}/${meta.normal_id}.bam_metric.summary.tsv" : ''

    if (tumor_arg && ! purple_arg) {
        error 'TEAL requires PURPLE inputs when analysing tumor data'
    }

    if (! tumor_arg && ! reference_arg) {
        error 'TEAL requires at least tumor or normal data for analyses'
    }

    """
    teal \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.teal.TealPipelineApp \\
        ${args} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${tumor_arg} \\
        ${tumor_bam_arg} \\
        -cobalt_dir ${cobalt_dir} \\
        ${purple_arg} \\
        ${reference_wgs_metrics_arg} \\
        ${tumor_wgs_metrics_arg} \\
        -ref_genome_version ${genome_ver} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir teal/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teal: \$(teal -version | sed -n '/Teal version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p teal/

    ${ (meta.tumor_id != null) ? "gzip <<< '' > teal/${meta.tumor_id}.teal.breakend.tsv.gz" : '' }
    ${ (meta.tumor_id != null) ? "touch teal/${meta.tumor_id}.teal.tellength.tsv" : '' }
    ${ (meta.normal_id != null) ? "touch teal/${meta.normal_id}.teal.tellength.tsv" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
