process SAGE_APPEND {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/sage:3.2.5--0'

    input:
    tuple val(meta), path(vcf), path(bam), path(bai)
    path genome_fasta
    path genome_fai
    path genome_dict

    output:
    tuple val(meta), path('*.append.vcf.gz'), emit: vcf
    path 'versions.yml'                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} com.hartwig.hmftools.sage.append.SageAppendApplication \\
            ${args} \\
            -input_vcf ${vcf} \\
            -reference ${meta.tumor_wts_id} \\
            -reference_bam ${bam} \\
            -ref_genome ${genome_fasta} \\
            -threads ${task.cpus} \\
            -out ./${meta.wgs_id}.sage.append.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*Sage version: //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.wgs_id}.sage.append.vcf.gz"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
