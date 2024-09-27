process PEACH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/hmftools-bam-tools:1.2.1--hdfd78af_0' :
//        'biocontainers/hmftools-bam-tools:1.2.1--hdfd78af_0' }"

    container 'quay.io/local/hmftools-peach'

    input:
    tuple val(meta), path(sage_germline_vcf), path(sage_germline_tbi)
    path(haplotypes)
    path(haplotype_functions)
    path(drugs)

    output:
    tuple val(meta), path('peach/'), emit: peach_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p peach/

    peach \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -vcf_file ${sage_germline_vcf} \\
        -sample_name ${meta.sample_id} \\
        -haplotypes_file ${haplotypes} \\
        -function_file ${haplotype_functions} \\
        -drugs_file ${drugs} \\
        -log_level DEBUG \\
        -output_dir peach/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peach: \$(peach -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch peach/${meta.sample_id}.peach.events.tsv
    touch peach/${meta.sample_id}.peach.gene.events.tsv
    touch peach/${meta.sample_id}.peach.haplotypes.all.tsv
    touch peach/${meta.sample_id}.peach.haplotypes.best.tsv
    touch peach/${meta.sample_id}.peach.qc.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
