nextflow.enable.types = true

process PURPLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-purple:4.4--hdfd78af_0' :
        'biocontainers/hmftools-purple:4.4--hdfd78af_0' }"

    input:
    tuple(meta: Map, amber_dir: Path, cobalt_dir: Path, esvee_dir: Path?, pave_somatic_dir: Path?, pave_germline_dir: Path?, redux_tumor_tsvs: List<Path>?)
    genome_fasta: Path
    genome_ver: String
    genome_fai: Path
    genome_dict: Path
    gc_profile: Path
    sage_known_hotspots_somatic: Path
    sage_known_hotspots_germline: Path?
    driver_gene_panel: Path
    ensembl_data_resources: Path
    germline_amp_del_freq: Path?
    target_regions_bed: Path?

    stage:
    stageAs redux_tumor_tsvs, 'redux_tumor_tsvs/*'

    topic:
    tuple(meta, file('purple/')) >> 'purple_dir'
    tuple(meta, 'purple', files('.command.*')) >> 'command_files'
    file('versions.yml') >> 'versions'

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
    def redux_tumor_dir_arg = redux_tumor_tsvs ? '-redux_tumor_dir redux_tumor_tsvs/' : ''

    def sage_known_hotspots_germline_arg = sage_known_hotspots_germline ? "-germline_hotspots ${sage_known_hotspots_germline}" : ''
    def germline_amp_del_freq_file_arg = germline_amp_del_freq ? "-germline_amp_del_freq_file ${germline_amp_del_freq}" : ''

    def target_regions_bed_arg = target_regions_bed ? "-target_regions_bed ${target_regions_bed}" : ''

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
        ${target_regions_bed_arg} \\
        ${germline_amp_del_freq_file_arg} \\
        -gc_profile ${gc_profile} \\
        -circos \$(which circos) \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir purple/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: \$(purple -version | sed -n '/^Purple version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        r-dplyr: \$(Rscript -e 'packageVersion("dplyr") |> as.character() |> writeLines()')
        r-ggplot2: \$(Rscript -e 'packageVersion("ggplot2") |> as.character() |> writeLines()')
        circos: \$(circos -version | sed -n '/^circos/ { s/^.* v //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p purple/ purple/plot/

    touch purple/plot/${meta.tumor_id}.circos.png

    touch purple/${meta.tumor_id}.purple.cnv.gene.tsv
    touch purple/${meta.tumor_id}.purple.cnv.somatic.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.germline.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.somatic.tsv
    gzip <<< '' > purple/${meta.tumor_id}.purple.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.germline.vcf.gz.tbi
    touch purple/${meta.tumor_id}.purple.purity.tsv
    touch purple/${meta.tumor_id}.purple.qc
    gzip <<< '' > purple/${meta.tumor_id}.purple.somatic.vcf.gz
    touch purple/${meta.tumor_id}.purple.somatic.vcf.gz.tbi
    gzip <<< '' >purple/${meta.tumor_id}.purple.sv.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.germline.vcf.gz.tbi
    gzip <<< '' > purple/${meta.tumor_id}.purple.sv.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
