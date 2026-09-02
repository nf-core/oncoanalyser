nextflow.enable.types = true

process NEO_ANNOTATE_FUSIONS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-isofox:2.0.1--hdfd78af_0' :
        'biocontainers/hmftools-isofox:2.0.1--hdfd78af_0' }"

    input:
    tuple(meta: Record, neo_finder_dir: Path, aln: Path, idx: Path)
    read_length: Integer
    genome_fasta: Path
    genome_ver: String
    genome_fai: Path
    ensembl_data_resources: Path

    topic:
    tuple(meta, file('*isf.neoepitope.tsv'))                 >> 'neo_annotated_fusions_tsv'
    tuple(meta, 'neo_annotate_fusions', files('.command.*')) >> 'command_files'
    file('versions.yml')                                     >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    isofox \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam_file ${aln} \\
        -functions NEO_EPITOPES \\
        -read_length ${read_length} \\
        -neo_dir ${neo_finder_dir} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox -version | sed -n '/^Isofox version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.isf.neoepitope.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
