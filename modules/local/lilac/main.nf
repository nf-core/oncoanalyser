process LILAC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-lilac:1.6--hdfd78af_1' :
        'biocontainers/hmftools-lilac:1.6--hdfd78af_1' }"

    input:
    tuple val(meta), path(normal_dna_bam), path(normal_dna_bai), path(tumor_dna_bam), path(tumor_dna_bai), path(tumor_rna_bam), path(tumor_rna_bai), path(purple_dir)
    path genome_fasta
    path genome_fai
    val genome_ver
    path lilac_resources, stageAs: 'lilac_resources'

    output:
    tuple val(meta), path('lilac/'), emit: lilac_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def sample_name = getSampleName(meta, tumor_dna_bam, normal_dna_bam)

    def normal_bam_arg = normal_dna_bam ? "-reference_bam ${normal_dna_bam}" : ''
    def tumor_dna_bam_arg = tumor_dna_bam ? "-tumor_bam ${tumor_dna_bam}" : ''
    def tumor_rna_bam_arg = tumor_rna_bam ? "-rna_bam ${tumor_rna_bam}" : ''

    def purple_dir_arg = purple_dir ? "-purple_dir ${purple_dir}" : ''

    """
    lilac \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${sample_name} \\
        ${normal_bam_arg} \\
        ${tumor_dna_bam_arg} \\
        ${tumor_rna_bam_arg} \\
        ${purple_dir_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -resource_dir ${lilac_resources} \\
        -threads ${task.cpus} \\
        -output_dir lilac/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lilac: \$(lilac -version | sed -n '/^Lilac version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p lilac/
    touch lilac/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

def getSampleName(meta, tumor_bam, normal_bam) {
    if (tumor_bam) {
        return meta.tumor_id
    } else if (normal_bam) {
        return meta.normal_id
    } else {
        Sys.exit(1)
    }
}
