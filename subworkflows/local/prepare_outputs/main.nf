//
// Prepare results to be published by entry workflow output block
//

// NOTE(SW): this approach is used so that linkage does not need to be maintained until the `output` block and absent outputs are handled implicitly and to avoid have verbose / repetitive emit blocks

workflow PREPARE_OUTPUTS_WGTS {
    main:
    ch_results = channel.empty()
        .mix(
            channel.topic('amber_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('bamtools_metrics_dir').flatMap { meta, d ->         return get_dir_filepaths(meta, d, "bamtools/${meta.sample_id}") },
            channel.topic('chord_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('cider_results').flatMap { meta, fps ->              return fps.collect { d -> ["${meta.key}/cider/${d.name}", d] } },
            channel.topic('cobalt_dir').flatMap { meta, d ->                   return get_dir_filepaths(meta, d) },
            channel.topic('cuppa_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('esvee_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('gatk4_markduplicates_bai').map { meta, d ->         return ["${meta.key}/alignments/${meta.sample_id}/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bam').map { meta, d ->         return ["${meta.key}/alignments/${meta.sample_id}/${d.name}", d] },
            channel.topic('isofox_dir').flatMap { meta, d ->                   return get_dir_filepaths(meta, d) },
            channel.topic('lilac_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('linx_germline_annotation_dir').flatMap { meta, d -> return get_dir_filepaths(meta, d, 'linx/germline_annotations') },
            channel.topic('linx_somatic_annotation_dir').flatMap { meta, d ->  return get_dir_filepaths(meta, d, 'linx/somatic_annotations') },
            channel.topic('linx_visualiser_plots').flatMap { meta, d ->        return get_dir_filepaths(meta, d, 'linx/somatic_plots') },
            channel.topic('linxreport_html').map { meta, d ->                  return ["${meta.key}/linx/${d.name}", d] },
            channel.topic('multiqc_report').map { d ->                         return [d.name, d] },
            channel.topic('neo_annotated_fusions_tsv').map { meta, d ->        return ["${meta.key}/neo/annotated_fusions/${d.name}", d] },
            channel.topic('neo_finder_dir').flatMap { meta, d ->               return get_dir_filepaths(meta, d, 'neo/finder') },
            channel.topic('neo_scorer_dir').flatMap { meta, d ->               return get_dir_filepaths(meta, d, 'neo/scorer') },
            channel.topic('orange_json').map { meta, d ->                      return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('orange_pdf').map { meta, d ->                       return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('pave_germline_dir').flatMap { meta, d ->            return get_dir_filepaths(meta, d, 'pave/germline') },
            channel.topic('pave_somatic_dir').flatMap { meta, d ->             return get_dir_filepaths(meta, d, 'pave/somatic') },
            channel.topic('peach_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('purple_dir').flatMap { meta, d ->                   return get_dir_filepaths(meta, d) },
            channel.topic('qsee_dir').flatMap { meta, d ->                     return get_dir_filepaths(meta, d) },
            channel.topic('redux_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d, "alignments/${meta.sample_id}") },
            channel.topic('sage_append_dir').flatMap { meta, d ->              return get_dir_filepaths(meta, d, "sage_append/${meta.output_file_id}") },
            channel.topic('sage_germline_dir').flatMap { meta, d ->            return get_dir_filepaths(meta, d, 'sage/germline') },
            channel.topic('sage_somatic_dir').flatMap { meta, d ->             return get_dir_filepaths(meta, d, 'sage/somatic') },
            channel.topic('sage_visualiser_dir').flatMap { meta, d ->          return get_dir_filepaths(meta, d, 'sage/visualiser') },
            channel.topic('sigs_dir').flatMap { meta, d ->                     return get_dir_filepaths(meta, d) },
            channel.topic('teal_prep_normal_bam').flatMap { meta, fps ->       return fps.collect { d -> ["${meta.key}/teal/${d.name}", d] } },
            channel.topic('teal_prep_tumor_bam').flatMap { meta, fps ->        return fps.collect { d -> ["${meta.key}/teal/${d.name}", d] } },
            channel.topic('teal_tsvs').flatMap { meta, e ->                    def fps = e instanceof List ? e : [e]; fps.collect { d -> ["${meta.key}/teal/${d.name}", d] } },
            channel.topic('virusbreakend_tsv').map { meta, d ->                return ["${meta.key}/virusbreakend/${d.name}", d] },
            channel.topic('virusbreakend_vcf').map { meta, d ->                return ["${meta.key}/virusbreakend/${d.name}", d] },
            channel.topic('virusinterpreter_dir').flatMap { meta, d ->         return get_dir_filepaths(meta, d) },

            channel.topic('write_reference_data').map { d -> return ["reference_data/${workflow.manifest.version}/", d] },

            channel.topic('command_files').flatMap { f -> get_command_log_filepath(f) }
        )
        .flatMap { meta, d -> return d instanceof List ? d.collect { e -> [meta, e] } : [[meta, d]] }

    emit:
    results = ch_results
}


workflow PREPARE_OUTPUTS_TARGETED {
    main:
    ch_results = channel.empty()
        .mix(
            channel.topic('amber_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('bamtools_metrics_dir').flatMap { meta, d ->         return get_dir_filepaths(meta, d, "bamtools/${meta.sample_id}") },
            channel.topic('cider_results').flatMap { meta, fps ->              return fps.collect { d -> ["${meta.key}/cider/${d.name}", d] } },
            channel.topic('cobalt_dir').flatMap { meta, d ->                   return get_dir_filepaths(meta, d) },
            channel.topic('esvee_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('gatk4_markduplicates_bai').map { meta, d ->         return ["${meta.key}/alignments/${meta.sample_id}/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bam').map { meta, d ->         return ["${meta.key}/alignments/${meta.sample_id}/${d.name}", d] },
            channel.topic('isofox_dir').flatMap { meta, d ->                   return get_dir_filepaths(meta, d) },
            channel.topic('lilac_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('linx_germline_annotation_dir').flatMap { meta, d -> return get_dir_filepaths(meta, d, 'linx/germline_annotations') },
            channel.topic('linx_somatic_annotation_dir').flatMap { meta, d ->  return get_dir_filepaths(meta, d, 'linx/somatic_annotations') },
            channel.topic('linx_visualiser_plots').flatMap { meta, d ->        return get_dir_filepaths(meta, d, 'linx/somatic_plots') },
            channel.topic('linxreport_html').map { meta, d ->                  return ["${meta.key}/linx/${d.name}", d] },
            channel.topic('multiqc_report').map { d ->                         return [d.name, d] },
            channel.topic('orange_json').map { meta, d ->                      return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('orange_pdf').map { meta, d ->                       return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('pave_germline_dir').flatMap { meta, d ->            return get_dir_filepaths(meta, d, 'pave/germline') },
            channel.topic('pave_somatic_dir').flatMap { meta, d ->             return get_dir_filepaths(meta, d, 'pave/somatic') },
            channel.topic('peach_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            channel.topic('purple_dir').flatMap { meta, d ->                   return get_dir_filepaths(meta, d) },
            channel.topic('qsee_dir').flatMap { meta, d ->                     return get_dir_filepaths(meta, d) },
            channel.topic('redux_dir').flatMap { meta, d ->                    return get_dir_filepaths(meta, d, "alignments/${meta.sample_id}") },
            channel.topic('sage_append_dir').flatMap { meta, d ->              return get_dir_filepaths(meta, d, "sage_append/${meta.output_file_id}") },
            channel.topic('sage_germline_dir').flatMap { meta, d ->            return get_dir_filepaths(meta, d, 'sage/germline') },
            channel.topic('sage_somatic_dir').flatMap { meta, d ->             return get_dir_filepaths(meta, d, 'sage/somatic') },
            channel.topic('sage_visualiser_dir').flatMap { meta, d ->          return get_dir_filepaths(meta, d, 'sage/visualiser') },

            channel.topic('write_reference_data').map { d -> return ["reference_data/${workflow.manifest.version}/", d] },

            channel.topic('command_files').flatMap { f -> get_command_log_filepath(f) }
        )
        .flatMap { meta, d -> return d instanceof List ? d.collect { e -> [meta, e] } : [[meta, d]] }

    emit:
    results = ch_results
}


workflow PREPARE_OUTPUTS_PURITY_ESTIMATE {
    main:
    ch_results = channel.empty()
        .mix(
            channel.topic('amber_dir').flatMap { meta, d ->            return get_dir_filepaths(meta, d) },
            channel.topic('cobalt_dir').flatMap { meta, d ->           return get_dir_filepaths(meta, d) },
            channel.topic('redux_dir').flatMap { meta, d ->            return get_dir_filepaths(meta, d, "alignments/${meta.sample_id}") },
            channel.topic('sage_append_dir').flatMap { meta, d ->      return get_dir_filepaths(meta, d, "sage_append/${meta.output_file_id}") },
            channel.topic('wisp_dir').flatMap { meta, d ->             return get_dir_filepaths(meta, d) },

            channel.topic('write_reference_data').map { d -> return ["reference_data/${workflow.manifest.version}/", d] },

            channel.topic('command_files').flatMap { f -> get_command_log_filepath(f) }
        )
        .flatMap { meta, d -> return d instanceof List ? d.collect { e -> [meta, e] } : [[meta, d]] }

    emit:
    results = ch_results
}


workflow PREPARE_OUTPUTS_PANEL_RESOURCE_CREATION {
    main:
    ch_results = channel.empty()
        .mix(
            channel.topic('amber_dir').flatMap { meta, d ->               return get_dir_filepaths(meta, d) },
            channel.topic('cobalt_dir').flatMap { meta, d ->              return get_dir_filepaths(meta, d) },
            channel.topic('gatk4_markduplicates_bai').map { meta, d ->    return ["${meta.key}/alignments/${meta.sample_id}/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bam').map { meta, d ->    return ["${meta.key}/alignments/${meta.sample_id}/${d.name}", d] },
            channel.topic('isofox_dir').flatMap { meta, d ->              return get_dir_filepaths(meta, d) },
            channel.topic('redux_dir').flatMap { meta, d ->               return get_dir_filepaths(meta, d, "alignments/${meta.sample_id}") },
            channel.topic('sage_germline_dir').flatMap { meta, d ->       return get_dir_filepaths(meta, d, 'sage/germline') },
            channel.topic('sage_somatic_dir').flatMap { meta, d ->        return get_dir_filepaths(meta, d, 'sage/somatic') },

            channel.topic('cobalt_normalisation_tsv').map { d ->          return ["panel_resources/${d.name}", d] },
            channel.topic('isofox_normalisation_csv').map { d ->          return ["panel_resources/${d.name}", d] },
            channel.topic('pave_pon_panel_creation_artefacts').map { d -> return ["panel_resources/${d.name}", d] },

            channel.topic('write_reference_data').map { d -> return ["reference_data/${workflow.manifest.version}/", d] },

            channel.topic('command_files').flatMap { f -> get_command_log_filepath(f) }
        )
        .flatMap { meta, d -> return d instanceof List ? d.collect { e -> [meta, e] } : [[meta, d]] }

    emit:
    results = ch_results
}


workflow PREPARE_OUTPUTS_PREPARE_REFERENCE {
    main:
    ch_results = channel.empty()
        .mix(
            channel.topic('write_reference_data').map { d -> return ["reference_data/${workflow.manifest.version}/", d] },

            channel.topic('command_files').flatMap { f -> get_command_log_filepath(f) }
        )

    emit:
    results = ch_results
}

// NOTE(LN): Nextflow bug as of NXF_VER=26.04.6
//
// When using google cloud executor (and potentially other cloud executors):
// - Process emits output as directory (e.g. amber_dir) -> directory publish FAIL -> BUG
// - Process emits output as files or glob (e.g. teal_tsvs) -> individual files publish SUCCESS
//
// To publish directories, the workaround is to recursively extract the file paths within that directory
//
def get_dir_filepaths(meta, d, target_dir=null) {

    def dir = target_dir ?: d.name

    def filepaths = []
    d.eachFileRecurse(groovy.io.FileType.FILES) { f ->
        def rel_path = d.relativize(f).toString()
        filepaths << ["${meta.key}/${dir}/${rel_path}", f]
    }
    return filepaths
}


def get_command_log_filepath(data) {

    def other_logs = ['gatk4_bwa_index_image', 'gridss_index', 'bwa_index', 'bwamem2_index', 'samtools_dict', 'samtools_faidx', 'star_genomegenerate', 'multiqc']
    def panel_logs = ['cobalt_panel_normalisation', 'pave_pon_panel_creation', 'isofox_panel_normalisation']

    def (meta, name, fps_all) = data

    def fps = fps_all.findAll { f -> f.name.matches(/.*\.command\.(sh|out|err|log|run)/) }

    if (name == 'extracttarball') {
        return fps.collect { d -> ["logs/other/${name}.${meta.id}${d.name}", d] }
    } else if (other_logs.contains(name)) {
        return fps.collect { d -> ["logs/other/${name}${d.name}", d] }
    } else if (panel_logs.contains(name)) {
        return fps.collect { d -> ["logs/panel_resources/${name}${d.name}", d] }
    } else {
        return fps.collect { d -> ["logs/${meta.key}/${name}.${meta.id}${d.name}", d] }
    }

}
