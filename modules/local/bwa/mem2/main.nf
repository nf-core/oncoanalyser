process BWA_MEM2 {
    tag "${meta.subject_id}__${meta.sample_id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--hd03093a_5' :
        'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    // TODO(SW): The following resourse files are needed from gs://hmf-public/HMFtools-Resources/ref_genome/37:
    //    + Homo_sapiens.GRCh37.GATK.illumina.fasta.bwt.2bit.64
    //    + Homo_sapiens.GRCh37.GATK.illumina.fasta.0123
    //    Similarly for ref genome 38.
    path genome_bwa_index

    output:
    tuple val(meta), path('*.bam'), emit: bam
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // TODO(MC): Double check this with Charles.
    def read_group_tag = "@RG\\tID:${meta.read_group}\\tSM:${meta.sample_id}"

    // NOTE(MC): Hardcoding bwa-mem2 version since the CLI does not have a --version flag.
    """
    ln -fs \$(find -L ${genome_bwa_index} -type f) ./

    bwa-mem2 mem \\
    -Y \\
    -R '${read_group_tag}' \\
    -t ${task.cpus} \\
    ${genome_fasta} \\
    ${reads_fwd} \\
    ${reads_rev} | \\
    \\
    sambamba view \\
      --sam-input \\
      --format bam \\
      --compression-level 0 \\
      --nthreads ${task.cpus} \\
      /dev/stdin | \\
    \\
    sambamba sort \\
      --nthreads ${task.cpus} \\
      --out ${meta.split}.${meta.sample_id}.${meta.read_group}.bam \\
      /dev/stdin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: 2.2.1
        sambamba: \$(sambamba --version 2>&1 | egrep '^sambamba' | head -n 1 | awk '{ print \$NF }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.split}.${meta.sample_id}.${meta.read_group}.bam

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}