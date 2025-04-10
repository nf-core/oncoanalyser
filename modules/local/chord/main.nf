process CHORD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-chord:2.1.0--hdfd78af_0' :
        'biocontainers/hmftools-chord:2.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(smlv_vcf), path(sv_vcf)
    path genome_fasta
    path genome_fai
    path genome_dict

    output:
    tuple val(meta), path('chord/'), emit: chord_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    """
    ## NOTE(LN): The CHORD jar runs an embedded R script using 'com.hartwig.hmftools.common.utils.r.RExecutor' which requires absolute
    ## paths. Relative paths don't work because RExecutor executes from a tmp dir, and not the working dir of this nextflow process

    mkdir -p chord/

    chord \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -snv_indel_vcf_file \$(realpath ${smlv_vcf}) \\
        -sv_vcf_file \$(realpath ${sv_vcf}) \\
        -output_dir \$(realpath chord/) \\
        -ref_genome ${genome_fasta} \\
        -log_level DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chord: \$(chord -version | sed -n '/^CHORD version/ { s/^.* //p }')
    END_VERSIONS

    """

    stub:
    """
    mkdir -p chord/
    touch chord/${meta.sample_id}.chord.mutation_contexts.tsv
    touch chord/${meta.sample_id}.chord.prediction.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
