nextflow.enable.types = true

process NEO_FINDER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-neo:1.3--hdfd78af_0' :
        'biocontainers/hmftools-neo:1.3--hdfd78af_0' }"

    input:
    tuple(meta: Record, purple_dir: Path, linx_annotation_dir: Path)
    genome_fasta: Path
    genome_ver: String
    genome_fai: Path
    ensembl_data_resources: Path

    topic:
    tuple(meta, file('neo_finder/'))               >> 'neo_finder_dir'
    tuple(meta, 'neo_finder', files('.command.*')) >> 'command_files'
    file('versions.yml')                           >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    mkdir -p neo_finder/

    neo \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -linx_dir ${linx_annotation_dir} \\
        -somatic_vcf ${purple_dir}/${meta.sample_id}.purple.somatic.vcf.gz \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        ${log_level_arg} \\
        -output_dir neo_finder/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neo: \$(neo -version | sed -n '/^Neo version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p neo_finder/

    touch neo_finder/.stub

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
