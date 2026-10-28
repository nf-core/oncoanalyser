process REDUX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-redux:2.0.5--hdfd78af_0' :
        'biocontainers/hmftools-redux:2.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(alns), path(idxs)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path unmap_regions
    path msi_jitter_sites
    path msi_model_coefficients
    path msi_model_error_rates
    val sequencing_platform
    val targeted_mode
    val generate_tsvs_only
    val umi_enable
    val umi_duplex_delim

    output:
    tuple val(meta), path("redux_${meta.sample_id}/"), topic: redux_dir
    tuple val(meta), val('redux'), path('.command.*'), topic: command_files
    path 'versions.yml'                              , topic: versions

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

    def bqr_use_all_regions_arg = ''

    def bqr_jitter_msi_only_arg = ''
    def aln_symlink_ext = ''
    def idx_symlink_ext = ''

    def msi_model_coefficients_arg = ''
    def msi_model_error_rates_arg = ''

    if (umi_enable) {
        umi_enable_arg = '-umi_enabled'
    } else {
        form_consensus_arg = '-form_consensus'
    }

    if (umi_duplex_delim) {
        umi_duplex_arg = '-umi_duplex'
        umi_duplex_delim_arg = "-umi_duplex_delim ${umi_duplex_delim}"
    }

    def umi_args = [umi_enable_arg, umi_duplex_arg, umi_duplex_delim_arg]
        .findAll { s -> s != '' }
        .join(' ')

    if (sequencing_platform.toLowerCase() == 'ultima') {
        form_consensus_arg = ''
        skip_duplicate_marking_arg = '-skip_duplicate_marking'
    }

    if (targeted_mode) {
        bqr_use_all_regions_arg = '-bqr_use_all_regions'
    }

    if (generate_tsvs_only) {
        def bad_inputs = [alns].flatten().size() != 1 || [idxs].flatten().size() != 1
        if (bad_inputs) { error 'Got too many alignments / indexes for `generate_tsvs_only`' }

        aln_symlink_ext = alns.extension
        idx_symlink_ext = idxs.extension

        bqr_jitter_msi_only_arg = '-bqr_jitter_msi_only'
    }

    if (meta.sample_type == 'tumor') {

        if (msi_model_coefficients) {
            msi_model_coefficients_arg = "-msi_model_coefficients ${msi_model_coefficients}"
        }

        if (msi_model_error_rates) {
            msi_model_error_rates_arg  = "-msi_model_error_rates ${msi_model_error_rates}"
        }

    }

    """
    mkdir -p redux_${meta.sample_id}/

    redux \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_bam ${alns.join(',')} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ref_genome_msi_file ${msi_jitter_sites} \\
        -unmap_regions ${unmap_regions} \\
        -bamtool \$(which samtools) \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -bqr_write_plot \\
        ${msi_model_coefficients_arg} \\
        ${msi_model_error_rates_arg} \\
        ${form_consensus_arg} \\
        ${umi_args} \\
        ${skip_duplicate_marking_arg} \\
        ${bqr_use_all_regions_arg} \\
        ${bqr_jitter_msi_only_arg} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_bam redux_${meta.sample_id}/${meta.sample_id}.redux.bam \\
        -output_dir redux_${meta.sample_id}/

    if [[ -n '${bqr_jitter_msi_only_arg}' ]]; then
        ln -sf \$(realpath ${alns}) redux_${meta.sample_id}/${meta.sample_id}.redux.${aln_symlink_ext}
        ln -sf \$(realpath ${idxs}) redux_${meta.sample_id}/${meta.sample_id}.redux.${aln_symlink_ext}.${idx_symlink_ext}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        redux: \$(redux -version | sed -n '/^Redux version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        r-dplyr: \$(Rscript -e 'packageVersion("dplyr") |> as.character() |> writeLines()')
        r-ggplot2: \$(Rscript -e 'packageVersion("ggplot2") |> as.character() |> writeLines()')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p redux_${meta.sample_id}/

    touch redux_${meta.sample_id}/${meta.sample_id}.redux.bam
    touch redux_${meta.sample_id}/${meta.sample_id}.redux.bam.bai

    touch redux_${meta.sample_id}/${meta.sample_id}.redux.bqr.tsv
    touch redux_${meta.sample_id}/${meta.sample_id}.redux.bqr.png
    touch redux_${meta.sample_id}/${meta.sample_id}.redux.duplicate_freq.tsv
    touch redux_${meta.sample_id}/${meta.sample_id}.redux.jitter_params.tsv
    touch redux_${meta.sample_id}/${meta.sample_id}.redux.msi_prediction.tsv
    gzip <<< '' > redux_${meta.sample_id}/${meta.sample_id}.redux.ms_table.tsv.gz
    gzip <<< '' > redux_${meta.sample_id}/${meta.sample_id}.redux.repeat.tsv.gz

    if [[ "${umi_enable}" == "true" ]]; then
        touch redux_${meta.sample_id}/${meta.sample_id}.umi_coord_freq.tsv
        touch redux_${meta.sample_id}/${meta.sample_id}.umi_edit_distance.tsv
        touch redux_${meta.sample_id}/${meta.sample_id}.umi_nucleotide_freq.tsv
    fi;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
