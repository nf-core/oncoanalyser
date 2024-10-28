process CHORD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-chord:2.1.0_beta--r43hdfd78af_0' :
        'biocontainers/hmftools-chord:2.1.0_beta--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(smlv_vcf), path(sv_vcf)
    val genome_ver

    output:
    tuple val(meta), path('chord/'), emit: chord_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p purple/
    ln -sf \$(realpath ${smlv_vcf}) purple/
    ln -sf \$(realpath ${sv_vcf}) purple/

    mkdir -p chord/

    ## NOTE(LN): Use pwd so that absolute path can be specified to -purple_dir and -output_dir
    ## Relative paths don't work because the RExecutor class from hmf-common executes from a tmp dir, and not the working dir of this
    ## nextflow process
    working_dir=\$PWD

    chord \\
        -Xmx${Math.round(task.memory.bytes * 0.75)} \\
        com.hartwig.hmftools.chord.ChordRunner \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -ref_genome_version ${genome_ver} \\
        -purple_dir "\${working_dir}/purple/" \\
        -output_dir "\${working_dir}/chord/" \\
        -log_level DEBUG

    ## NOTE(LN): Chord expects the R packages `CHORD` and `mutSigExtractor` to be installed as R packages
    ## Otherwise, arg '-chord_tool_dir' would need to be provided, whereby this dir contains the subdirs CHORD/ and mutSigExtractor/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chord: \$(chord -version | awk '{ print \$NF }')
    END_VERSIONS

    """

    stub:
    """
    mkdir -p chord/
    touch chord/${meta.sample_id}_chord_signatures.txt
    touch chord/${meta.sample_id}_chord_prediction.txt

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
