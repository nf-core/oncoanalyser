process ESVEE_PREP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.0_beta--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-esvee:1.0_beta--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path genome_fasta
    val genome_ver
    path sv_blocklist
    path known_fusions

    output:
    tuple val(meta), path("sv_prep/")                                                , emit: sv_prep_dir
    tuple val(meta), path("sv_prep/${meta.normal_id}.esvee.prep.bam")                , emit: normal_prep_bam
    tuple val(meta), path("sv_prep/${meta.tumor_id}.esvee.prep.bam")                 , emit: tumor_prep_bam
    tuple val(meta), path("sv_prep/${meta.tumor_id}.esvee.prep.junctions.tsv")       , emit: junctions_tsv
    tuple val(meta), path("sv_prep/${meta.tumor_id}.esvee.prep.fragment_lengths.tsv"), emit: fragment_lengths_tsv
    path 'versions.yml'                                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def sample_ids = [meta.tumor_id]
    def bam_files = [tumor_bam.toString()]

    if(meta.normal_id != null){
        sample_ids.add(meta.normal_id)
        bam_files.add(normal_bam.toString())
    }

    def sample_ids_string = String.join(",", sample_ids)
    def bam_files_string = String.join(",", bam_files)

    """
    mkdir -p sv_prep/

    SAMBAMBA_PATH=\$(which sambamba)

    esvee com.hartwig.hmftools.esvee.prep.SvPrepApplication \\
        -Xmx${Math.round(task.memory.bytes * 0.75)} \\
        ${args} \\
        -sample "${sample_ids_string}" \\
        -bam_files "${bam_files_string}" \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -blacklist_bed ${sv_blocklist} \\
        -known_fusion_bed ${known_fusions} \\
        -bamtool \$SAMBAMBA_PATH \\
        -write_types "JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST" \\
        -output_dir sv_prep/ \\
        -threads ${task.cpus} \\
        -log_level DEBUG \\

    # NOTE(LN): For tumor only mode, make empty output null.esvee.prep.bam for the reference sample so that nextflow doesn't complain about
    # this missing file when emitting output
    ${ (meta.normal_id == null) ? "touch sv_prep/${meta.normal_id}.esvee.prep.bam" : "" }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sv_prep/

    touch "sv_prep/${meta.normal_id}.esvee.prep.bam"
    touch "sv_prep/${meta.normal_id}.esvee.prep.bam.bai"
    touch "sv_prep/${meta.tumor_id}.esvee.prep.bam"
    touch "sv_prep/${meta.tumor_id}.esvee.prep.bam.bai"
    touch "sv_prep/${meta.tumor_id}.esvee.prep.fragment_lengths.tsv"
    touch "sv_prep/${meta.tumor_id}.esvee.prep.junctions.tsv"

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
