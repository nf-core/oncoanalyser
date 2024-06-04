process WISP {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/wisp:1.2_beta.2--0'

    input:
    tuple val(meta), path(sage_append_dir), path('sample_amber_dir'), path(cobalt_dir), path('primary_amber_dir'), path(primary_purple_dir)
    path genome_fasta
    path genome_fai

    output:
    path 'wisp/'       , emit: wisp_dir
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def somatic_dir_arg = sage_append_dir ? "-somatic_dir somatic_dir__prepared/" : ''
    def amber_dir_arg = sample_amber_dir && primary_amber_dir ? "-amber_dir amber_dir__prepared/" : ''
    def cobalt_dir_arg = cobalt_dir ? "-cobalt_dir ${cobalt_dir}" : ''
    def purple_dir_arg = primary_purple_dir ? "-purple_dir ${primary_purple_dir}" : ''

    def purity_methods_list = []
    if (somatic_dir_arg) purity_methods_list.add('SOMATIC_VARIANT')
    if (amber_dir_arg) purity_methods_list.add('AMBER_LOH')
    if (cobalt_dir_arg) purity_methods_list.add('COPY_NUMBER')

    def purity_methods_arg = purity_methods_list ? "-purity_methods '${purity_methods_list.join(';')}'" : ''

    """
    # Prepare SAGE append directory (-somatic_dir)
    if [[ -n "${somatic_dir_arg}" ]]; then
      mkdir -p somatic_dir__prepared/
      ln -s ../${sage_append_dir}/${meta.primary_id}.sage.append.vcf.gz somatic_dir__prepared/${meta.primary_id}.purple.somatic.ctdna.vcf.gz;
      ln -s ../${sage_append_dir}/${meta.primary_id}.sage.append.vcf.gz.tbi somatic_dir__prepared/${meta.primary_id}.purple.somatic.ctdna.vcf.gz.tbi;
      ln -s ../${sage_append_dir}/${meta.sample_id}.sage.bqr.tsv somatic_dir__prepared/;
      ln -s ../${sage_append_dir}/${meta.sample_id}.frag_lengths.tsv.gz somatic_dir__prepared/;
    fi

    # Prepare AMBER directory (-amber_dir)
    if [[ -n "${amber_dir_arg}" ]]; then
      mkdir -p amber_dir__prepared/;
      for fp in ${primary_amber_dir}/${meta.primary_id}.*; do ln -fs ../\${fp} amber_dir__prepared/; done;
      for fp in ${sample_amber_dir}/${meta.sample_id}.*; do fn=\${fp##*/}; ln -fs ../\${fp} amber_dir__prepared/; done;
    fi;

    # Run WISP
    mkdir -p wisp/

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp /opt/wisp/wisp.jar com.hartwig.hmftools.wisp.purity.PurityEstimator \\
            ${args} \\
            -patient_id ${meta.subject_id} \\
            -tumor_id ${meta.primary_id} \\
            -samples ${meta.sample_id} \\
            \\
            ${purity_methods_arg} \\
            \\
            ${somatic_dir_arg} \\
            ${amber_dir_arg} \\
            ${cobalt_dir_arg} \\
            ${purple_dir_arg} \\
            \\
            -ref_genome ${genome_fasta} \\
            \\
            -write_types ALL \\
            -log_debug \\
            \\
            -output_dir wisp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisp: \$(java -jar /opt/wisp/wisp.jar -version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p wisp/
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.cn_plot_calcs.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.cn_segments.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.somatic_peak.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.somatic_variants.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.summary.tsv
    touch wisp/${meta.sample_id}.cn_gc_ratio_fit.png
    touch wisp/${meta.sample_id}.somatic_vaf.png

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
