process ESVEE_DEPTH_ANNOTATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.0_beta--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-esvee:1.0_beta--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(raw_vcf)
    path genome_fasta
    val genome_ver

    output:
    tuple val(meta), path("depth_annotation/")                                           , emit: depth_annotation_dir
    tuple val(meta), path("depth_annotation/${meta.tumor_id}.esvee.ref_depth.vcf.gz")    , emit: ref_depth_vcf
    path 'versions.yml'                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def sample_ids = [meta.tumor_id]
    def bam_files = [tumor_bam.toString()]

    if(meta.normal_id != null){
        sample_ids.add(meta.normal_id)
        bam_files.add(normal_bam.toString())
    }

    def sample_ids_string = String.join(",", sample_ids)
    def bam_files_string = String.join(",", bam_files)

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
    touch depth_annotation/${meta.tumor_id}.esvee.ref_depth.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
