process ESVEE_PREP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.4--hdfd78af_0' :
        'docker.io/scwatts/hmftools-esvee:1.0_beta--hdfd78af_0--1' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path genome_fasta
    val genome_ver
    path sv_blocklist
    path known_fusions

    output:
    tuple val(meta), path('sv_prep/'), emit: sv_prep_dir
    path 'versions.yml'              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def sample_ids_string = String.join(",", meta.tumor_id, meta.normal_id)
    def bam_files_string = String.join(",", tumor_bam.toString(), normal_bam.toString())

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
