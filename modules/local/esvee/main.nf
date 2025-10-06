process ESVEE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.1.2--hdfd78af_0' :
        'biocontainers/hmftools-esvee:1.1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_img
    val genome_ver
    path pon_breakends
    path pon_breakpoints
    path decoy_sequences_image
    path known_fusions
    path repeatmasker_annotations
    path unmap_regions
    path target_region_bed

    output:
    tuple val(meta), path('esvee/')                                                                                                    , emit: esvee_dir
    tuple val(meta), path("esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz"), path("esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf
    tuple val(meta), path("esvee/${meta.tumor_id}.esvee.somatic.vcf.gz"),    path("esvee/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi")   , emit: somatic_vcf
    tuple val(meta), path("esvee/${meta.tumor_id}.esvee.germline.vcf.gz"),   path("esvee/${meta.tumor_id}.esvee.germline.vcf.gz.tbi")  , emit: germline_vcf, optional: true
    path 'versions.yml'                                                                                                                , emit: versions
    path '.command.*'                                                                                                                  , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_arg = meta.normal_id ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = meta.normal_id ? "-reference_bam ${normal_bam}" : ''

    def target_region_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''

    """
    mkdir -p esvee/

    esvee \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -esvee_prep_dir esvee/ \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -known_hotspot_file ${known_fusions} \\
        -pon_sgl_file ${pon_breakends} \\
        -pon_sv_file ${pon_breakpoints} \\
        -repeat_mask_file ${repeatmasker_annotations} \\
        -unmap_regions ${unmap_regions} \\
        ${target_region_bed_arg} \\
        -bamtool \$(which sambamba) \\
        -write_types 'PREP_JUNCTION;PREP_BAM;FRAGMENT_LENGTH_DIST;JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENT;BREAKEND;VCF' \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir esvee/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(java -jar \${ESVEE_JAR} -version | sed 's/^.*Esvee version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p esvee/

    touch esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz
    touch esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi
    touch esvee/${meta.tumor_id}.esvee.somatic.vcf.gz
    touch esvee/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi
    touch esvee/${meta.tumor_id}.esvee.germline.vcf.gz
    touch esvee/${meta.tumor_id}.esvee.germline.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(echo "1.0-beta")
    END_VERSIONS
    """
}
