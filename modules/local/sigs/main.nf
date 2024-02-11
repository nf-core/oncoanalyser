process SIGS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sigs:1.2.1--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-sigs:1.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(smlv_vcf)
    path signatures

    output:
    tuple val(meta), path('sigs/'), emit: sigs_dir
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p sigs/

    sigs \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -sample ${meta.sample_id} \\
        -somatic_vcf_file ${smlv_vcf} \\
        -signatures_file ${signatures} \\
        -output_dir sigs/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigs: \$(sigs -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sigs/
    touch sigs/placeholder
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
