process ESVEE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.1.2--hdfd78af_0' :
        'biocontainers/hmftools-esvee:1.1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_prep_bam), path(tumor_prep_bai), path(normal_prep_bam), path(normal_prep_bai)
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_img
    val genome_ver
    path pon_breakends
    path pon_breakpoints
    path known_fusions
    path repeatmasker_annotations
    path unmap_regions

    output:
    tuple val(meta), path("esvee/"), emit: esvee_dir
    tuple val(meta), path("esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz"), path("esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf
    tuple val(meta), path("esvee/${meta.tumor_id}.esvee.somatic.vcf.gz"),    path("esvee/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi"),    emit: somatic_vcf
    tuple val(meta), path("esvee/${meta.tumor_id}.esvee.germline.vcf.gz"),   path("esvee/${meta.tumor_id}.esvee.germline.vcf.gz.tbi"),   emit: germline_vcf, optional: true
    path 'versions.yml', emit: versions
    path '.command.*'  , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.normal_id ? "-reference ${meta.normal_id}" : ""
    def tumor_bam_arg = "-tumor_bam ${tumor_prep_bam}"

    def reference_bam_arg = meta.normal_id ? "-reference_bam ${normal_prep_bam}" : ""

    def sample_ids = [meta.tumor_id]

    if (meta.normal_id) {
        sample_ids.add(meta.normal_id)
    }


    """
    mkdir -p esvee/

    esvee \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        ${reference_arg} \\
        ${tumor_bam_arg} \\
        ${reference_bam_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -known_hotspot_file ${known_fusions} \\
        -pon_sgl_file ${pon_breakends} \\
        -pon_sv_file ${pon_breakpoints} \\
        -repeat_mask_file ${repeatmasker_annotations} \\
        -unmap_regions ${unmap_regions} \\
        -bamtool \$(which sambamba) \\
        -write_types "PREP_JUNCTION;PREP_BAM;FRAGMENT_LENGTH_DIST;JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENT;BREAKEND;VCF" \\
        -output_dir esvee/ \\
        -esvee_prep_dir esvee/ \\
        -threads ${task.cpus} \\
        -log_level ${params.module_log_level}

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
