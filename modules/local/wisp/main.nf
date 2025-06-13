process WISP {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/wisp:1.2_beta.2--0'

    input:
    tuple val(meta),
        path(primary_purple_dir),
        path('primary_amber_dir'),
        path('sample_amber_dir'),
        path(cobalt_dir),
        path(sage_append_dir)
    path genome_fasta
    path genome_fai
    val is_targeted_mode

    output:
    path 'wisp/'       , emit: wisp_dir
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def purity_estimate_mode = Utils.getEnumFromString(params.purity_estimate_mode, Constants.RunMode)

    def purity_methods
    def amber_dir_arg
    def cobalt_dir_arg
    def gc_ratio_min_arg

    if(is_targeted_mode) {
        purity_methods      = "'SOMATIC_VARIANT'"
        amber_dir_arg       = ""
        cobalt_dir_arg      = ""
        gc_ratio_min_arg    = "-gc_ratio_min 0.4"
    } else {
        purity_methods      = "'SOMATIC_VARIANT;AMBER_LOH;COPY_NUMBER'"
        amber_dir_arg       = "-amber_dir amber_dir__prepared/"
        cobalt_dir_arg      = "-cobalt_dir ${cobalt_dir}"
        gc_ratio_min_arg    = ""
    }

    """
    # Rename SAGE append files to the names expected by WISP
    mkdir -p somatic_dir__prepared/
    ln -sf ../${sage_append_dir}/${meta.primary_id}.sage.append.vcf.gz              somatic_dir__prepared/${meta.primary_id}.purple.somatic.ctdna.vcf.gz;
    ln -sf ../${sage_append_dir}/${meta.primary_id}.sage.append.vcf.gz.tbi          somatic_dir__prepared/${meta.primary_id}.purple.somatic.ctdna.vcf.gz.tbi;
    ln -sf ../${sage_append_dir}/${meta.primary_id}.sage.append.frag_lengths.tsv.gz somatic_dir__prepared/${meta.primary_id}.frag_lengths.tsv.gz;

    # Put AMBER outputs from all samples into the same dir
    if [[ -n "${amber_dir_arg}" ]]; then
        mkdir -p amber_dir__prepared/;
        for fp in ${primary_amber_dir}/*.amber.*; do ln -sf ../\$fp amber_dir__prepared/; done
        for fp in ${sample_amber_dir}/*.amber.*;  do ln -sf ../\$fp amber_dir__prepared/; done
    fi;

    # Run WISP
    mkdir -p wisp/

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp /opt/wisp/wisp.jar com.hartwig.hmftools.wisp.purity.PurityEstimator \\
        ${args} \\
        -patient_id ${meta.subject_id} \\
        -tumor_id ${meta.primary_id} \\
        -samples ${meta.longitudinal_id} \\
        -ref_genome ${genome_fasta} \\
        -purity_methods ${purity_methods} \\
        -somatic_dir somatic_dir__prepared/ \\
        -purple_dir ${primary_purple_dir} \\
        ${amber_dir_arg} \\
        ${cobalt_dir_arg} \\
        ${gc_ratio_min_arg} \\
        -output_dir wisp/ \\
        -write_types ALL \\
        -log_level ${params.module_log_level}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisp: \$(java -jar /opt/wisp/wisp.jar -version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p wisp/
    touch wisp/${meta.subject_id}_${meta.longitudinal_id}.wisp.cn_plot_calcs.tsv
    touch wisp/${meta.subject_id}_${meta.longitudinal_id}.wisp.cn_segments.tsv
    touch wisp/${meta.subject_id}_${meta.longitudinal_id}.wisp.somatic_peak.tsv
    touch wisp/${meta.subject_id}_${meta.longitudinal_id}.wisp.somatic_variants.tsv
    touch wisp/${meta.subject_id}_${meta.longitudinal_id}.wisp.summary.tsv
    touch wisp/${meta.longitudinal_id}.cn_gc_ratio_fit.png
    touch wisp/${meta.longitudinal_id}.somatic_vaf.png

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
