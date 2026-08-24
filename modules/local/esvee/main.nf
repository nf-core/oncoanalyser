nextflow.enable.types = true

include { NormalReferenceMeta } from '../../../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/records'

process ESVEE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:2.0.1--hdfd78af_0' :
        'biocontainers/hmftools-esvee:2.0.1--hdfd78af_0' }"

    input:
    tuple(meta: NormalReferenceMeta, tumor_aln: Path, tumor_bai: Path, normal_aln: Path?, normal_bai: Path?)
    genome_fasta: Path
    genome_ver: String
    genome_fai: Path
    genome_dict: Path
    genome_img: Path
    pon_breakends: Path
    pon_breakpoints: Path
    decoy_sequences_image: Path?
    known_fusions: Path
    repeatmasker_annotations: Path
    unmap_regions: Path
    target_regions_bed: Path?
    sequencing_platform: String

    topic:
    tuple(meta, file('esvee/'))               >> 'esvee_dir'
    tuple(meta, 'esvee', files('.command.*')) >> 'command_files'
    file('versions.yml')                      >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_arg = meta.normal_id != null ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = meta.normal_id != null ? "-reference_bam ${normal_aln}" : ''

    def decoy_genome_arg = decoy_sequences_image ? "-decoy_genome ${decoy_sequences_image}" : ''

    def target_regions_bed_arg = target_regions_bed ? "-target_regions_bed ${target_regions_bed}" : ''

    """
    mkdir -p esvee/

    esvee \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_aln} \\
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
        ${decoy_genome_arg} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        ${target_regions_bed_arg} \\
        -bamtool \$(which sambamba) \\
        -write_types 'PREP_JUNCTION;PREP_BAM;FRAGMENT_LENGTH_DIST;JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENT;BREAKEND;VCF' \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir esvee/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed -n '/^.*Esvee version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        sambamba: \$(sambamba --version 2>&1 | sed -n '/^sambamba / { s/^.* //p }' | head -n1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p esvee/

    gzip <<< '' > esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz
    touch esvee/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi
    gzip <<< '' > esvee/${meta.tumor_id}.esvee.somatic.vcf.gz
    touch esvee/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi
    gzip <<< '' > esvee/${meta.tumor_id}.esvee.germline.vcf.gz
    touch esvee/${meta.tumor_id}.esvee.germline.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
