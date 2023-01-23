process PURPLE {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/purple:3.7.1--0'

    input:
    tuple val(meta), path(amber), path(cobalt), path(sv_hard_vcf), path(sv_hard_vcf_index), path(sv_soft_vcf), path(sv_soft_vcf_index), path(smlv_tumor_vcf), path(smlv_normal_vcf)
    path genome_fasta
    path genome_fai
    path genome_dict
    val genome_ver
    path gc_profile
    path sage_known_hotspots_somatic
    path sage_known_hotspots_germline
    path driver_gene_panel
    path ensembl_data_resources
    path germline_del

    output:
    tuple val(meta), path('purple/'), emit: purple_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def smlv_tumor_vcf_fp = smlv_tumor_vcf ?: ''
    def smlv_normal_vcf_fp = smlv_normal_vcf ?: ''
    def germline_del_arg = germline_del ? "-germline_del_freq_file ${germline_del}" : ''

    """
    # For provided smlv VCFs, filter records that do not contain the required FORMAT/AD field and
    # get argument for PURPLE
    get_smlv_arg() {
        fp=\${1}
        fn=\${fp##*/}
        if [[ "\${fp}" != '' ]]; then
            fp_out="prepared__\${2}__\${fn}"
            bcftools filter -Oz -e 'FORMAT/AD[*]="."' "\${fp}" > \${fp_out}
            echo "-\${2} \${fp_out}"
        fi
    }
    smlv_tumor_vcf_arg=\$(get_smlv_arg "${smlv_tumor_vcf_fp}" somatic_vcf)
    smlv_normal_vcf_arg=\$(get_smlv_arg "${smlv_normal_vcf_fp}" germline_vcf)

    # Run PURPLE
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -tumor ${meta.tumor_id} \\
            -reference ${meta.normal_id} \\
            -sv_recovery_vcf ${sv_soft_vcf} \\
            -structural_vcf ${sv_hard_vcf} \\
            \${smlv_tumor_vcf_arg} \\
            \${smlv_normal_vcf_arg} \\
            -amber ${amber} \\
            -cobalt ${cobalt} \\
            -output_dir purple/ \\
            -gc_profile ${gc_profile} \\
            -run_drivers \\
            -driver_gene_panel ${driver_gene_panel} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -somatic_hotspots ${sage_known_hotspots_somatic} \\
            -germline_hotspots ${sage_known_hotspots_germline} \\
            ${germline_del_arg} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -threads ${task.cpus} \\
            -circos ${task.ext.circosPath}

    # PURPLE can fail silently, check that at least the PURPLE SV VCF is created
    if [[ ! -s "purple/${meta.tumor_id}.purple.sv.vcf.gz" ]]; then
        exit 1;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: \$(java -jar ${task.ext.jarPath} -version | sed 's/.*Purple version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir purple/
    touch purple/${meta.tumor_id}.purple.cnv.gene.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.germline.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.somatic.tsv
    touch purple/${meta.tumor_id}.purple.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.purity.tsv
    touch purple/${meta.tumor_id}.purple.qc
    touch purple/${meta.tumor_id}.purple.somatic.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
