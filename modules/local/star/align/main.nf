nextflow.enable.types = true

process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.3a--0' :
        'biocontainers/star:2.7.3a--0' }"

    input:
    tuple(meta: Map, rg_lines: List<String>, reads_fwds: List<Path>, reads_revs: List<Path>)
    genome_star_index: Path

    topic:
    tuple(meta, file('*bam'))                      >> 'star_align_bam'
    tuple(meta, file('*Log.final.out'))            >> 'star_align_qc_log'
    tuple(meta, 'star_align', files('.command.*')) >> 'command_files'
    file('versions.yml')                           >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reads_fwds_str = reads_fwds.join(',')
    def reads_revs_str = reads_revs.join(',')

    def rg_lines_str = rg_lines.join(' , ')

    """
    STAR \\
        ${args} \\
        --readFilesIn ${reads_fwds_str} ${reads_revs_str} \\
        --outSAMattrRGline ${rg_lines_str} \\
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
    touch Log.final.out

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
