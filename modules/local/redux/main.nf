process REDUX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-redux:1.2.2--hdfd78af_0' :
        'biocontainers/hmftools-redux:1.2.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bams), path(bais)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path unmap_regions
    path msi_jitter_sites
    val sequencing_type
    val sequencing_type_ultima
    val umi_enable
    val umi_duplex_delim
    val targeted_mode

    output:
    tuple val(meta), path('*.redux.bam'),
                     path('*.redux.bam.bai')             , emit: bam

    tuple val(meta), path('*.redux.bqr.tsv'),
                     path('*.redux.duplicate_freq.tsv'),
                     path('*.redux.jitter_params.tsv'),
                     path('*.redux.ms_table.tsv.gz')     , emit: tsv

    tuple val(meta), path('*.redux.bqr.png')             , emit: bqr_plot
    path 'versions.yml'                                  , emit: versions
    path '.command.*'                                    , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def form_consensus_arg = ''
    def umi_enable_arg = ''
    def umi_duplex_arg = ''
    def umi_duplex_delim_arg = ''
    def skip_duplicate_marking_arg = ''

    if(umi_enable) {
        umi_enable_arg = '-umi_enabled'
    } else {
        form_consensus_arg = '-form_consensus'
    }

    if(umi_duplex_delim) {
        umi_duplex_arg = '-umi_duplex'
        umi_duplex_delim_arg = "-umi_duplex_delim ${umi_duplex_delim}"
    }

    def umi_args = [umi_enable_arg, umi_duplex_arg, umi_duplex_delim_arg]
        .findAll { it != '' }
        .join(' ')

    if(sequencing_type_ultima) {
        form_consensus_arg = ''
        skip_duplicate_marking_arg = '-skip_duplicate_marking'
    }

    def bqr_use_all_regions_arg = targeted_mode ? '-bqr_use_all_regions' : ''

    """
    redux \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_bam ${bams.join(',')} \\
        -output_bam ./${meta.sample_id}.redux.bam \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ref_genome_msi_file ${msi_jitter_sites} \\
        -unmap_regions ${unmap_regions} \\
        -bamtool \$(which samtools) \\
        -sequencing_type ${sequencing_type} \\
        -bqr_write_plot \\
        ${form_consensus_arg} \\
        ${umi_args} \\
        ${skip_duplicate_marking_arg} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        redux: \$(redux -version | sed -n '/^Redux version/ { s/^.* //p }')
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.redux.bam
    touch ${meta.sample_id}.redux.bam.bai
    touch ${meta.sample_id}.redux.bqr.tsv
    touch ${meta.sample_id}.redux.duplicate_freq.tsv
    touch ${meta.sample_id}.redux.jitter_params.tsv
    touch ${meta.sample_id}.redux.ms_table.tsv.gz
    touch ${meta.sample_id}.redux.repeat.tsv.gz

    if [[ -n "${umi_enable}" ]]; then
        touch ${meta.sample_id}.umi_coord_freq.tsv
        touch ${meta.sample_id}.umi_edit_distance.tsv
        touch ${meta.sample_id}.umi_nucleotide_freq.tsv
    fi;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
