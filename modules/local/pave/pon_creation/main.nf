process PAVE_PON_PANEL_CREATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.7_beta--hdfd78af_1' :
        'biocontainers/hmftools-pave:1.7--hdfd78af_0' }"

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
       echo SampleId;
       find ${sage_vcf} | sed 's#.sage.somatic.vcf.gz##';
    ) > sample_ids.txt

    java -cp /usr/local/share/hmftools-pave-1.7-0/pave.jar \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.pave.resources.PonBuilder \\
            ${args} \\
            -sample_id_file sample_ids.txt \\
            -vcf_path '*.sage.somatic.vcf.gz' \\
            -ref_genome_version ${genome_ver} \\
            -output_dir ./

    mv somatic_pon_*.tsv pave.somatic_artefacts.${genome_ver}.tsv

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

