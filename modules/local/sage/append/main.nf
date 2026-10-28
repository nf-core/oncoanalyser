process SAGE_APPEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:5.0.2--hdfd78af_0' :
        'biocontainers/hmftools-sage:5.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(alns), path(idxs), path(redux_tsvs)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    val sequencing_platform
    val targeted_mode

    output:
    tuple val(meta), path("sage_append_${meta.output_file_id}/"), topic: sage_append_dir
    tuple val(meta), val('sage_append'), path('.command.*')     , topic: command_files
    path 'versions.yml'                                         , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def skip_args = ! redux_tsvs ? '-skip_bqr -skip_msi_jitter' : ''
    def high_depth_mode_arg = targeted_mode ? '-high_depth_mode' : ''

    """
    mkdir -p sage_append_${meta.output_file_id}/

    sage \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.sage.append.SageAppendApplication \\
        ${args} \\
        -input_vcf ${vcf} \\
        -max_read_depth 100000 \\
        -reference ${meta.reference_ids.join(',')} \\
        -reference_bam ${alns.join(',')} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -write_frag_lengths \\
        ${high_depth_mode_arg} \\
        ${skip_args} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_vcf sage_append_${meta.output_file_id}/${meta.output_file_id}.sage.append.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed -n '/^Sage version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sage_append_${meta.output_file_id}/

    gzip <<< '' > sage_append_${meta.output_file_id}/${meta.output_file_id}.frag_lengths.tsv.gz
    gzip <<< '' > sage_append_${meta.output_file_id}/${meta.output_file_id}.sage.append.vcf.gz
    touch sage_append_${meta.output_file_id}/${meta.output_file_id}.sage.append.vcf.gz.tbi
    touch sage_append_${meta.output_file_id}/${meta.output_file_id}_query.sage.bqr.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
