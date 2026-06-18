process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.3a--0' :
        'biocontainers/star:2.7.3a--0' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_star_index

    output:
    tuple val(meta), path('*bam'), emit: bam
    path 'versions.yml'          , emit: versions
    path '.command.*'            , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // Build the read group line. Extra/override SAM fields may be supplied via the 'read_group' samplesheet
    // info field (e.g. 'ID=..|PL=..'); ID and SM default to the sample/library/lane identifiers. STAR expects
    // space-separated 'KEY:VALUE' attributes with no '@RG' prefix.
    def rg_fields = meta.read_group_fields ?: [:]
    def rg_id = rg_fields.ID ?: meta.read_group
    def rg_sm = rg_fields.SM ?: meta.sample_id
    def rg_extra = rg_fields.findAll { k, v -> k != 'ID' && k != 'SM' }.collect { k, v -> "${k}:${v}" }
    def read_group_line = (["ID:${rg_id}", "SM:${rg_sm}"] + rg_extra).join(' ')

    """
    STAR \\
        ${args} \\
        --readFilesIn ${reads_fwd} ${reads_rev} \\
        --genomeDir ${genome_star_index} \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --alignSplicedMateMapLmin 35 \\
        --alignSplicedMateMapLminOverLmate 0.33 \\
        --chimJunctionOverhangMin 10 \\
        --chimOutType WithinBAM SoftClip \\
        --chimScoreDropMax 70 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreMin 1 \\
        --chimScoreSeparation 1 \\
        --chimSegmentMin 10 \\
        --chimSegmentReadGapMax 3 \\
        --limitOutSJcollapsed 3000000 \\
        --outBAMcompression 0 \\
        --outFilterMatchNmin 35 \\
        --outFilterMatchNminOverLread 0.33 \\
        --outFilterMismatchNmax 3 \\
        --outFilterMultimapNmax 10 \\
        --outFilterScoreMinOverLread 0.33 \\
        --outSAMattributes All \\
        --outSAMattrRGline ${read_group_line} \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --peOverlapNbasesMin 10 \\
        --runRNGseed 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    """
    touch Aligned.out.bam

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
