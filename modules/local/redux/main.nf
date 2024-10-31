process REDUX {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-redux:1.0_beta--hdfd78af_2' :
        'biocontainers/hmftools-redux:1.0_beta--hdfd78af_2' }"

    input:
    tuple val(meta), path(bams), path(bais)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path unmap_regions
    path msi_jitter_sites
    val umi_enable
    val umi_duplex_delim

    output:
    tuple val(meta), path('*.redux.bam'), path('*.redux.bam.bai'), emit: bam
    tuple val(meta), path('*.duplicate_freq.tsv')                , emit: dup_freq_tsv
    tuple val(meta), path('*.jitter_params.tsv')                 , emit: jitter_tsv
    tuple val(meta), path('*.ms_table.tsv.gz')                   , emit: ms_tsv
    tuple val(meta), path('*.repeat.tsv.gz')                     , emit: repeat_tsv
    path 'versions.yml'                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def form_consensus_arg = umi_enable ? '' : '-form_consensus'

    def umi_args_list = []
    if (umi_enable) umi_args_list.add('-umi_enabled')
    if (umi_duplex_delim) umi_args_list.add("-umi_duplex -umi_duplex_delim ${umi_duplex_delim}")
    def umi_args = umi_args_list ? umi_args_list.join(' ') : ''

    """
    redux \\
        -Xmx${Math.round(task.memory.bytes * 0.75)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_bam ${bams.join(',')} \\
        -output_dir ./ \\
        -output_bam ./${meta.sample_id}.redux.bam \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -unmap_regions ${unmap_regions} \\
        -ref_genome_msi_file ${msi_jitter_sites} \\
        -bamtool \$(which sambamba) \\
        ${form_consensus_arg} \\
        ${umi_args} \\
        -write_stats \\
        -threads ${task.cpus} \\
        -log_level DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        redux: \$(redux -version | awk '{ print \$NF }')
        sambamba: \$(sambamba --version 2>&1 | egrep '^sambamba' | head -n 1 | awk '{ print \$NF }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.redux.bam
    touch ${meta.sample_id}.redux.bam.bai
    touch ${meta.sample_id}.duplicate_freq.tsv
    touch ${meta.sample_id}.jitter_params.tsv
    touch ${meta.sample_id}.ms_table.tsv.gz
    touch ${meta.sample_id}.repeat.tsv.gz

    if [[ -n "${umi_enable}" ]]; then
        touch ${meta.sample_id}.umi_coord_freq.tsv
        touch ${meta.sample_id}.umi_edit_distance.tsv
        touch ${meta.sample_id}.umi_nucleotide_freq.tsv
    fi;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
