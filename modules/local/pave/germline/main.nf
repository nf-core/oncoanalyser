// NOTE(SW): use of tumor sample name here is consistent with Pipeline5
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.33/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveGermline.java#L36-L41
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.33/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L31-L43

process PAVE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.7.1--hdfd78af_0' :
        'biocontainers/hmftools-pave:1.7.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_tbi)
    path genome_fasta
    val genome_ver
    path genome_fai
    path sage_blocklist_regions
    path sage_blocklist_sites
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def gnomad_args
    if (genome_ver.toString() == '37') {
        gnomad_args = "-gnomad_freq_file ${gnomad_resource}"
    } else if (genome_ver.toString() == '38') {
        gnomad_args = "-gnomad_freq_dir ${gnomad_resource}"
    } else {
        error "got bad genome version: ${genome_ver}"
    }

    """
    pave \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -vcf_file ${sage_vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -clinvar_vcf ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -mappability_bed ${segment_mappability} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -blacklist_bed ${sage_blocklist_regions} \\
        -blacklist_vcf ${sage_blocklist_sites} \\
        ${gnomad_args} \\
        -gnomad_no_filter \\
        -read_pass_only \\
        -threads ${task.cpus} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed -n '/^Pave version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.sage.pave_germline.vcf.gz{,.tbi}

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
