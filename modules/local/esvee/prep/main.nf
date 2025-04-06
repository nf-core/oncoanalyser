process ESVEE_PREP {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.0.3--hdfd78af_0' :
        'biocontainers/hmftools-esvee:1.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path genome_fasta
    val genome_ver
    path sv_blocklist
    path known_fusions

    output:
    tuple val(meta), path("prep/")                                                                 , emit: prep_dir
    tuple val(meta), path("prep/${meta.tumor_id}.*.bam"), path("prep/${meta.tumor_id}.*.bam.bai")  , emit: tumor_prep_bam
    tuple val(meta), path("prep/${meta.normal_id}.*.bam"), path("prep/${meta.normal_id}.*.bam.bai"), emit: normal_prep_bam, optional: true
    path 'versions.yml'                                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def sample_ids = [meta.tumor_id]
    def bam_files = [tumor_bam.toString()]

    if(meta.normal_id != null){
        sample_ids.add(meta.normal_id)
        bam_files.add(normal_bam.toString())
    }

    def sample_ids_string = String.join(',', sample_ids)
    def bam_files_string = String.join(',', bam_files)

    """
    mkdir -p prep/

    SAMBAMBA_PATH=\$(which sambamba)

    esvee com.hartwig.hmftools.esvee.prep.PrepApplication \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample "${sample_ids_string}" \\
        -bam_file "${bam_files_string}" \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -blacklist_bed ${sv_blocklist} \\
        -known_fusion_bed ${known_fusions} \\
        -bamtool \$SAMBAMBA_PATH \\
        -write_types 'JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST' \\
        -output_dir prep/ \\
        -threads ${task.cpus} \\
        -log_level DEBUG \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed -n '/^Esvee version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p prep/

    ${ (meta.normal_id != null) ? "touch prep/${meta.normal_id}.esvee.prep.bam" : '' }
    ${ (meta.normal_id != null) ? "touch prep/${meta.normal_id}.esvee.prep.bam.bai" : '' }
    touch "prep/${meta.tumor_id}.esvee.prep.bam"
    touch "prep/${meta.tumor_id}.esvee.prep.bam.bai"
    touch "prep/${meta.tumor_id}.esvee.prep.fragment_length.tsv"
    touch "prep/${meta.tumor_id}.esvee.prep.junction.tsv"

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
