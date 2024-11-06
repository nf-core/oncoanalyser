process PEACH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-peach:2.0.0--hdfd78af_1' :
        'biocontainers/hmftools-peach:2.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(germline_vcf)
    path haplotypes
    path haplotype_functions
    path drug_info

    output:
    tuple val(meta), path('peach/'), emit: peach_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    peach \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample_name ${meta.sample_id} \\
        -vcf_file ${germline_vcf} \\
        -haplotypes_file ${haplotypes} \\
        -function_file ${haplotype_functions} \\
        -drugs_file ${drug_info} \\
        -output_dir peach/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peach: \$(peach -version | sed -n '/Peach version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p peach/

    touch peach/${meta.sample_id}.peach.events.tsv
    touch peach/${meta.sample_id}.peach.gene.events.tsv
    touch peach/${meta.sample_id}.peach.haplotypes.all.tsv
    touch peach/${meta.sample_id}.peach.haplotypes.best.tsv
    touch peach/${meta.sample_id}.peach.qc.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
