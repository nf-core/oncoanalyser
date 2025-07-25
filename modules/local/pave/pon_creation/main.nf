process PAVE_PON_PANEL_CREATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.8--hdfd78af_1' :
        'biocontainers/hmftools-pave:1.8--hdfd78af_1' }"

    input:
    tuple path(sage_vcf), path(sage_tbi)
    val genome_ver

    output:
    path 'pave.somatic_artefacts.*.tsv', emit: pave_artefacts
    path 'versions.yml'                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    (
       echo SampleId
       basename -s .sage.somatic.vcf.gz -a *.sage.somatic.vcf.gz
    ) > sample_ids.txt

    pave \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.pave.pon_gen.PonBuilder \\
        ${args} \\
        -sample_id_file sample_ids.txt \\
        -vcf_path '*.sage.somatic.vcf.gz' \\
        -ref_genome_version ${genome_ver} \\
        -output_pon_file pave.somatic_artefacts.${genome_ver}.tsv \\
        -log_level ${params.module_log_level}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch pave.somatic_artefacts.${genome_ver}.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

