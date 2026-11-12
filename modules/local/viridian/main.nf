process VIRIDIAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // TODO: confirm once hmftools-viridian 1.0 is released on bioconda
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-viridian:1.0--hdfd78af_0' :
        'biocontainers/hmftools-viridian:1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_aln), path(tumor_idx), path(esvee_dir), path(bamtools_somatic_dir), path(purple_dir)
    path genome_fasta
    path genome_fai
    path genome_dict
    path viridian_resources, stageAs: 'viridian_resources'

    output:
    tuple val(meta), path('viridian/')                  , topic: viridian_dir
    tuple val(meta), val('viridian'), path('.command.*'), topic: command_files
    path 'versions.yml'                                 , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    mkdir -p viridian/

    viridian \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -tumor_bam ${tumor_aln} \\
        -esvee_unfiltered_vcf ${esvee_dir}/${meta.sample_id}.esvee.unfiltered.vcf.gz \\
        -purple_dir ${purple_dir} \\
        -bam_metrics_tumor_dir ${bamtools_somatic_dir} \\
        -ref_genome ${genome_fasta} \\
        -virus_genomes_fasta ${viridian_resources}/virus_genomes.fasta \\
        -virus_genomes_tsv ${viridian_resources}/virus_genomes.tsv \\
        -virus_groups_tsv ${viridian_resources}/virus_groups.tsv \\
        -kmer_filter \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir viridian/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viridian: \$(viridian -version | sed -n '/^Viridian version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p viridian/

    touch viridian/${meta.sample_id}.viridian.report.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
