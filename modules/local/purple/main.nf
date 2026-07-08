process PURPLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-purple:4.4--hdfd78af_0' :
        'biocontainers/hmftools-purple:4.4--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(amber_dir),
        path(cobalt_dir),
        path(esvee_dir),
        path(pave_somatic_dir, stageAs: 'pave_somatic'),
        path(pave_germline_dir, stageAs: 'pave_germline'),
        path(redux_tumor_tsvs, stageAs: "redux_tumor/*")
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path gc_profile
    path sage_known_hotspots_somatic
    path sage_known_hotspots_germline
    path driver_gene_panel
    path ensembl_data_resources
    path germline_amp_del_freq
    path target_region_bed

    output:
    tuple val(meta), path('purple/'), emit: purple_dir
    path 'versions.yml'             , emit: versions
    path '.command.*'               , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''

    def esvee_dir_arg = esvee_dir ? "-esvee_dir ${esvee_dir}" : ''
    def pave_somatic_dir_arg = pave_somatic_dir ? "-pave_somatic_dir ${pave_somatic_dir}" : ''
    def pave_germline_dir_arg = pave_germline_dir ? "-pave_germline_dir ${pave_germline_dir}" : ''
    def redux_tumor_dir_arg = redux_tumor_tsvs ? "-redux_tumor_dir redux_tumor/" : ''

    def sage_known_hotspots_germline_arg = sage_known_hotspots_germline ? "-germline_hotspots ${sage_known_hotspots_germline}" : ''
    def germline_amp_del_freq_file_arg = germline_amp_del_freq ? "-germline_amp_del_freq_file ${germline_amp_del_freq}" : ''

    def target_region_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''

    """
    purple \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        ${reference_arg} \\
        -amber ${amber_dir} \\
        -cobalt ${cobalt_dir} \\
        ${esvee_dir_arg} \\
        ${pave_somatic_dir_arg} \\
        ${pave_germline_dir_arg} \\
        ${redux_tumor_dir_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -somatic_hotspots ${sage_known_hotspots_somatic} \\
        ${sage_known_hotspots_germline_arg} \\
        ${target_region_bed_arg} \\
        ${germline_amp_del_freq_file_arg} \\
        -gc_profile ${gc_profile} \\
        -circos \$(which circos) \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir purple/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: \$(purple -version | sed -n '/^Purple version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir purple/

    touch purple/${meta.tumor_id}.purple.cnv.gene.tsv
    touch purple/${meta.tumor_id}.purple.cnv.somatic.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.germline.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.somatic.tsv
    touch purple/${meta.tumor_id}.purple.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.germline.vcf.gz.tbi
    touch purple/${meta.tumor_id}.purple.purity.tsv
    touch purple/${meta.tumor_id}.purple.qc
    touch purple/${meta.tumor_id}.purple.somatic.vcf.gz
    touch purple/${meta.tumor_id}.purple.somatic.vcf.gz.tbi
    touch purple/${meta.tumor_id}.purple.sv.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.germline.vcf.gz.tbi
    touch purple/${meta.tumor_id}.purple.sv.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
