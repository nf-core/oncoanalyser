process PAVE_PON_PANEL_CREATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.9--hdfd78af_0' :
        'biocontainers/hmftools-pave:1.9--hdfd78af_0' }"

    input:
    tuple path(sage_vcf), path(sage_tbi)
    val genome_ver

    output:
    path 'pave.somatic_artefacts.*.tsv'                               , topic: pave_pon_panel_creation_artefacts
    tuple val([:]), val('pave_pon_panel_creation'), path('.command.*'), topic: command_files
    path 'versions.yml'                                               , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

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
        ${log_level_arg} \\
        -output_pon_file pave.somatic_artefacts.${genome_ver}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed -n '/^Pave version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    touch pave.somatic_artefacts.${genome_ver}.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
