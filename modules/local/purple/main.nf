process PURPLE {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/purple:3.8.2--0'

    input:
    tuple val(meta), path(amber), path(cobalt), path(sv_tumor_vcf), path(sv_tumor_tbi), path(sv_tumor_unfiltered_vcf), path(sv_tumor_unfiltered_tbi), path(sv_normal_vcf), path(sv_normal_tbi), path(smlv_tumor_vcf), path(smlv_normal_vcf)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path gc_profile
    path sage_known_hotspots_somatic
    path sage_known_hotspots_germline
    path driver_gene_panel
    path ensembl_data_resources
    path germline_del
    path target_region_bed
    path target_region_ratios
    path target_region_msi_indels

    output:
    tuple val(meta), path('purple/'), emit: purple_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''

    def sv_tumor_vcf_arg = sv_tumor_vcf ? "-somatic_sv_vcf ${sv_tumor_vcf}" : ''
    def sv_normal_vcf_arg = sv_normal_vcf ? "-germline_sv_vcf ${sv_normal_vcf}" : ''

    def sv_tumor_recovery_vcf_arg = sv_tumor_unfiltered_vcf ? "-sv_recovery_vcf ${sv_tumor_unfiltered_vcf}" : ''

    def smlv_tumor_vcf_fp = smlv_tumor_vcf ?: ''
    def smlv_normal_vcf_fp = smlv_normal_vcf ?: ''

    def sage_known_hotspots_germline_arg = sage_known_hotspots_germline ? "-germline_hotspots ${sage_known_hotspots_germline}" : ''
    def germline_del_arg = germline_del ? "-germline_del_freq_file ${germline_del}" : ''

    def target_region_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''
    def target_region_ratios_arg = target_region_ratios ? "-target_regions_ratios ${target_region_ratios}" : ''
    def target_region_msi_indels_arg = target_region_msi_indels ? "-target_regions_msi_indels ${target_region_msi_indels}" : ''

    """
    # For provided smlv VCFs, filter records that do not contain the required FORMAT/AD field and
    # get argument for PURPLE
    get_smlv_arg() {
        fp=\${1}
        fn=\${fp##*/}
        if [[ "\${fp}" != '' ]]; then
            fp_out=prepared__\${2}__\${fn}
            bcftools view -e 'FORMAT/AD[*]="."' -o \${fp_out} \${fp}
            echo "-\${2} \${fp_out}"
        fi
    }
    smlv_tumor_vcf_arg=\$(get_smlv_arg "${smlv_tumor_vcf_fp}" somatic_vcf)
    smlv_normal_vcf_arg=\$(get_smlv_arg "${smlv_normal_vcf_fp}" germline_vcf)

    # Run PURPLE
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -tumor ${meta.tumor_id} \\
            ${reference_arg} \\
            -amber ${amber} \\
            -cobalt ${cobalt} \\
            ${sv_tumor_vcf_arg} \\
            ${sv_normal_vcf_arg} \\
            ${sv_tumor_recovery_vcf_arg} \\
            \${smlv_tumor_vcf_arg} \\
            \${smlv_normal_vcf_arg} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -somatic_hotspots ${sage_known_hotspots_somatic} \\
            ${sage_known_hotspots_germline_arg} \\
            ${target_region_bed_arg} \\
            ${target_region_ratios_arg} \\
            ${target_region_msi_indels_arg} \\
            ${germline_del_arg} \\
            -gc_profile ${gc_profile} \\
            -circos ${task.ext.circosPath} \\
            -threads ${task.cpus} \\
            -output_dir purple/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: \$(java -jar ${task.ext.jarPath} -version | sed 's/.*Purple version: //')
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
    touch purple/${meta.tumor_id}.purple.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.purity.tsv
    touch purple/${meta.tumor_id}.purple.qc
    touch purple/${meta.tumor_id}.purple.somatic.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
