process SAGE_APPEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.4--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sage:3.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf), path(bam), path(bai)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict

    output:
    tuple val(meta), path('*.append.vcf.gz'), emit: vcf
    path 'versions.yml'                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.sage.append.SageAppendApplication \\
        ${args} \\
        -input_vcf ${vcf} \\
        -reference ${meta.tumor_rna_id} \\
        -reference_bam ${bam} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -threads ${task.cpus} \\
        -output_vcf ${meta.dna_id}.sage.append.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.dna_id}.sage.append.vcf.gz"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
