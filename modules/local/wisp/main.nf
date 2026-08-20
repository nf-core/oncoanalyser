nextflow.enable.types = true

process WISP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-wisp:1.3.1--hdfd78af_0' :
        'biocontainers/hmftools-wisp:1.3.1--hdfd78af_0' }"

    input:
    tuple(meta: Map, primary_purple_dir: Path, primary_amber_dir: Path?, primary_normal_aln: Path?, longitudinal_redux_dir: Path, longitudinal_amber_dir: Path, longitudinal_cobalt_dir: Path, longitudinal_sage_append_dir: Path)
    genome_fasta: Path
    genome_fai: Path
    targeted_mode: Boolean

    stage:
    stageAs primary_purple_dir, 'purple_primary'
    stageAs primary_amber_dir, 'amber_primary'
    stageAs longitudinal_redux_dir, 'redux_longitudinal'
    stageAs longitudinal_amber_dir, 'amber_longitudinal'
    stageAs longitudinal_cobalt_dir, 'cobalt_longitudinal'
    stageAs longitudinal_sage_append_dir, 'sage_append_longitudinal'

    topic:
    tuple(meta, file('wisp/')) >> 'wisp_dir'
    tuple(meta, 'wisp', files('.command.*')) >> 'command_files'
    file('versions.yml') >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def amber_dir_arg
    def cobalt_dir_arg
    def gc_ratio_min_arg
    def write_types_arg

    def purity_methods = ['SOMATIC_VARIANT']

    if (targeted_mode) {

        amber_dir_arg = ''
        cobalt_dir_arg = ''
        gc_ratio_min_arg = '-gc_ratio_min 0.4'
        write_types_arg = "-write_types 'SOMATIC_DATA;SOMATIC_PLOT'"

    } else {

        if (primary_amber_dir && primary_normal_aln) {
            amber_dir_arg = '-amber_dir amber_dir__prepared/'
            purity_methods += 'AMBER_LOH'
        } else {
            amber_dir_arg = ''
        }

        cobalt_dir_arg = "-cobalt_dir ${longitudinal_cobalt_dir}"
        purity_methods += 'COPY_NUMBER'

        gc_ratio_min_arg = ''
        write_types_arg = '-write_types ALL'
    }

    def purity_methods_arg = "'${purity_methods.join(';')}'"

    """
    ## Put AMBER outputs from primary and longitudinal sample into the same dir
    if [[ -n "${amber_dir_arg}" ]]; then
        mkdir -p amber_dir__prepared/;
        for fp in ${primary_amber_dir}/*.amber.*; do ln -sf ../\$fp amber_dir__prepared/; done
        for fp in ${longitudinal_amber_dir}/*.amber.*; do ln -sf ../\$fp amber_dir__prepared/; done
    fi;

    mkdir -p wisp/

    wisp \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.wisp.purity.PurityEstimator \\
        ${args} \\
        -patient_id ${meta.patient_id} \\
        -tumor_id ${meta.primary_id} \\
        -samples ${meta.longitudinal_id} \\
        -purity_methods ${purity_methods_arg} \\
        -somatic_vcf ${longitudinal_sage_append_dir}/${meta.longitudinal_id}.sage.append.vcf.gz \\
        -purple_dir ${primary_purple_dir} \\
        -bqr_dir ${longitudinal_redux_dir} \\
        ${amber_dir_arg} \\
        ${cobalt_dir_arg} \\
        -ref_genome ${genome_fasta} \\
        ${gc_ratio_min_arg} \\
        ${write_types_arg} \\
        ${log_level_arg} \\
        -output_dir wisp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisp: \$(wisp -version | sed -n '/^Wisp version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        r-tidyverse: \$(Rscript -e 'packageVersion("tidyverse") |> as.character() |> writeLines()')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p wisp/

    touch wisp/${meta.patient_id}_${meta.longitudinal_id}.wisp.cn_plot_calcs.tsv
    touch wisp/${meta.patient_id}_${meta.longitudinal_id}.wisp.cn_segments.tsv
    touch wisp/${meta.patient_id}_${meta.longitudinal_id}.wisp.somatic_peak.tsv
    touch wisp/${meta.patient_id}_${meta.longitudinal_id}.wisp.somatic_variants.tsv
    touch wisp/${meta.patient_id}_${meta.longitudinal_id}.wisp.summary.tsv
    touch wisp/${meta.longitudinal_id}.cn_gc_ratio_fit.png
    touch wisp/${meta.longitudinal_id}.somatic_vaf.png

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
