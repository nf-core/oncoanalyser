process ESVEE_DEPTH_ANNOTATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.4--hdfd78af_0' :
        'docker.io/scwatts/hmftools-esvee:1.0_beta--hdfd78af_0--1' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai), path(raw_vcf)
    path genome_fasta
    val genome_ver

    output:
    tuple val(meta), path("depth_annotation/${meta.tumor_id}.esvee.ref_depth.vcf.gz")    , emit: ref_depth_vcf
    tuple val(meta), path("depth_annotation/${meta.tumor_id}.esvee.ref_depth.vcf.gz.tbi"), emit: ref_depth_tbi
    path 'versions.yml'                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def sample_ids_string = String.join(",", meta.tumor_id, meta.normal_id)
    def bam_files_string = String.join(",", tumor_bam.toString(), normal_bam.toString())

    """
    mkdir -p depth_annotation/

    esvee com.hartwig.hmftools.esvee.depth.DepthAnnotator \\
        -Xmx${Math.round(task.memory.bytes * 0.75)} \\
        ${args} \\
        -samples ${sample_ids_string} \\
        -bam_files ${bam_files_string} \\
        -input_vcf ${raw_vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -output_dir depth_annotation/ \\
        -threads ${task.cpus} \\
        -log_debug

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p depth_annotation/

    touch depth_annotation/${meta.tumor_id}.esvee.ref_depth.vcf.gz
    tough depth_annotation/${meta.tumor_id}.esvee.ref_depth.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
