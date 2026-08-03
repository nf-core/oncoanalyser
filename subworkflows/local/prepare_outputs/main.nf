//
// Prepare results to be published by entry workflow output block
//

// NOTE(SW): this approach is used so that linkage does not need to be maintained until the `output` block and absent outputs are handled implicitly and to avoid have verbose / repetitive emit blocks

workflow PREPARE_OUTPUTS_WGTS {
    main:
    ch_results = channel.empty()
        .mix(
            channel.topic('amber_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('bamtools_metrics_dir').map { meta, d ->         return ["${meta.key}/bamtools/${meta.sample_id}", d] },
            channel.topic('chord_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('cider_results').flatMap { meta, fps ->          return fps.collect { d -> ["${meta.key}/cider/${d.name}", d] } },
            channel.topic('cobalt_dir').map { meta, d ->                   return ["${meta.key}/${d.name}", d] },
            channel.topic('cuppa_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('esvee_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bai').map { meta, d ->     return ["${meta.key}/alignments/rna/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bam').map { meta, d ->     return ["${meta.key}/alignments/rna/${d.name}", d] },
            channel.topic('isofox_dir').map { meta, d ->                   return ["${meta.key}/${d.name}", d] },
            channel.topic('lilac_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('linx_germline_annotation_dir').map { meta, d -> return ["${meta.key}/linx/germline_annotations", d] },
            channel.topic('linx_somatic_annotation_dir').map { meta, d ->  return ["${meta.key}/linx/somatic_annotations", d] },
            channel.topic('linx_visualiser_plots').map { meta, d ->        return ["${meta.key}/linx/somatic_plots", d] },
            channel.topic('linxreport_html').map { meta, d ->              return ["${meta.key}/linx/${d.name}", d] },
            channel.topic('multiqc_report').map { d ->                     return [d.name, d] },
            channel.topic('neo_annotated_fusions_tsv').map { meta, d ->    return ["${meta.key}/neo/annotated_fusions/${d.name}", d] },
            channel.topic('neo_finder_dir').map { meta, d ->               return ["${meta.key}/neo/finder", d] },
            channel.topic('neo_scorer_dir').map { meta, d ->               return ["${meta.key}/neo/scorer", d] },
            channel.topic('orange_json').map { meta, d ->                  return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('orange_pdf').map { meta, d ->                   return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('pave_somatic_dir').map { meta, d ->             return ["${meta.key}/pave/somatic", d] },
            channel.topic('pave_germline_dir').map { meta, d ->            return ["${meta.key}/pave/germline", d] },
            channel.topic('peach_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('purple_dir').map { meta, d ->                   return ["${meta.key}/${d.name}", d] },
            channel.topic('qsee_dir').map { meta, d ->                     return ["${meta.key}/${d.name}", d] },
            channel.topic('redux_dir').flatMap { e -> def meta = e[0];     return e[1..-1].collect { d -> ["${meta.key}/alignments/dna/${d.name}", d] } },
            channel.topic('sage_append_dir').map { meta, d ->              return ["${meta.key}/sage_append/${meta.output_file_id}", d] },
            channel.topic('sage_germline_dir').map { meta, d ->            return ["${meta.key}/sage/germline", d] },
            channel.topic('sage_somatic_dir').map { meta, d ->             return ["${meta.key}/sage/somatic", d] },
            channel.topic('sage_visualiser_dir').map { meta, d ->          return ["${meta.key}/sage/visualiser", d] },
            channel.topic('sigs_dir').map { meta, d ->                     return ["${meta.key}/${d.name}", d] },
            channel.topic('teal_prep_normal_bam').flatMap { meta, fps ->   return fps.collect { d -> ["${meta.key}/teal/${d.name}", d] } },
            channel.topic('teal_prep_tumor_bam').flatMap { meta, fps ->    return fps.collect { d -> ["${meta.key}/teal/${d.name}", d] } },
            channel.topic('teal_tsvs').flatMap { meta, e ->                def fps = e instanceof List ? e : [e]; fps.collect { d -> ["${meta.key}/teal/${d.name}", d] } },
            channel.topic('virusbreakend_tsv').map { meta, d ->            return ["${meta.key}/virusbreakend/${d.name}", d] },
            channel.topic('virusbreakend_vcf').map { meta, d ->            return ["${meta.key}/virusbreakend/${d.name}", d] },
            channel.topic('virusinterpreter_dir').map { meta, d ->         return ["${meta.key}/${d.name}", d] },

            channel.topic('write_reference_data').map { d ->               return ["reference_data/${workflow.manifest.version}/", d] },

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
            channel.topic('amber_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('bamtools_metrics_dir').map { meta, d ->         return ["${meta.key}/bamtools/${meta.sample_id}", d] },
            channel.topic('cider_results').flatMap { meta, fps ->          return fps.collect { d -> ["${meta.key}/cider/${d.name}", d] } },
            channel.topic('cobalt_dir').map { meta, d ->                   return ["${meta.key}/${d.name}", d] },
            channel.topic('esvee_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bai').map { meta, d ->     return ["${meta.key}/alignments/rna/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bam').map { meta, d ->     return ["${meta.key}/alignments/rna/${d.name}", d] },
            channel.topic('isofox_dir').map { meta, d ->                   return ["${meta.key}/${d.name}", d] },
            channel.topic('lilac_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('linx_germline_annotation_dir').map { meta, d -> return ["${meta.key}/linx/germline_annotations", d] },
            channel.topic('linx_somatic_annotation_dir').map { meta, d ->  return ["${meta.key}/linx/somatic_annotations", d] },
            channel.topic('linx_visualiser_plots').map { meta, d ->        return ["${meta.key}/linx/somatic_plots", d] },
            channel.topic('linxreport_html').map { meta, d ->              return ["${meta.key}/linx/${d.name}", d] },
            channel.topic('multiqc_report').map { d ->                     return [d.name, d] },
            channel.topic('orange_json').map { meta, d ->                  return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('orange_pdf').map { meta, d ->                   return ["${meta.key}/orange/${d.name}", d] },
            channel.topic('pave_somatic_dir').map { meta, d ->             return ["${meta.key}/pave/somatic", d] },
            channel.topic('pave_germline_dir').map { meta, d ->            return ["${meta.key}/pave/germline", d] },
            channel.topic('peach_dir').map { meta, d ->                    return ["${meta.key}/${d.name}", d] },
            channel.topic('purple_dir').map { meta, d ->                   return ["${meta.key}/${d.name}", d] },
            channel.topic('redux_dir').flatMap { e -> def meta = e[0];     return e[1..-1].collect { d -> ["${meta.key}/alignments/dna/${d.name}", d] } },
            channel.topic('sage_append_dir').map { meta, d ->              return ["${meta.key}/sage_append/${meta.output_file_id}", d] },
            channel.topic('sage_germline_dir').map { meta, d ->            return ["${meta.key}/sage/germline", d] },
            channel.topic('sage_somatic_dir').map { meta, d ->             return ["${meta.key}/sage/somatic", d] },

            channel.topic('write_reference_data').map { d ->               return ["reference_data/${workflow.manifest.version}/", d] },

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
            channel.topic('amber_dir').map { meta, d ->                return ["${meta.key}/${d.name}", d] },
            channel.topic('cobalt_dir').map { meta, d ->               return ["${meta.key}/${d.name}", d] },
            channel.topic('redux_dir').flatMap { e -> def meta = e[0]; return e[1..-1].collect { d -> ["${meta.key}/alignments/dna/${d.name}", d] } },
            channel.topic('sage_append_dir').map { meta, d ->          return ["${meta.key}/sage_append/${meta.output_file_id}", d] },
            channel.topic('wisp_dir').map { meta, d ->                 return ["${meta.key}/${d.name}", d] },

            channel.topic('write_reference_data').map { d ->       return ["reference_data/${workflow.manifest.version}/", d] },

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
            channel.topic('amber_dir').map { meta, d ->                   return ["${meta.key}/${d.name}", d] },
            channel.topic('cobalt_dir').map { meta, d ->                  return ["${meta.key}/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bai').map { meta, d ->    return ["${meta.key}/alignments/rna/${d.name}", d] },
            channel.topic('gatk4_markduplicates_bam').map { meta, d ->    return ["${meta.key}/alignments/rna/${d.name}", d] },
            channel.topic('isofox_dir').map { meta, d ->                  return ["${meta.key}/${d.name}", d] },
            channel.topic('redux_dir').flatMap { e -> def meta = e[0];    return e[1..-1].collect { d -> ["${meta.key}/alignments/dna/${d.name}", d] } },
            channel.topic('sage_germline_dir').map { meta, d ->           return ["${meta.key}/sage/germline", d] },
            channel.topic('sage_somatic_dir').map { meta, d ->            return ["${meta.key}/sage/somatic", d] },

            channel.topic('cobalt_normalisation_tsv').map { d ->          return ["panel_resources/${d.name}", d] },
            channel.topic('isofox_normalisation_csv').map { d ->          return ["panel_resources/${d.name}", d] },
            channel.topic('pave_pon_panel_creation_artefacts').map { d -> return ["panel_resources/${d.name}", d] },

            channel.topic('write_reference_data').map { d ->              return ["reference_data/${workflow.manifest.version}/", d] },

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


def get_command_log_filepath(data) {

    def other_logs = ['gatk4_bwa_index_image', 'gridss_index', 'bwa_index', 'bwamem2_index', 'samtools_dict', 'samtools_faidx', 'star_genomegenerate', 'extracttarball', 'multiqc']
    def panel_logs = ['cobalt_panel_normalisation', 'pave_pon_panel_creation', 'isofox_panel_normalisation']

    def (meta, name, fps_all) = data

    def fps = fps_all.findAll { f -> f.name.matches(/.*\.command\.(sh|out|err|log|run)/) }

    if (other_logs.contains(name)) {
        return fps.collect { d -> ["logs/other/${name}${d.name}", d] }
    } else if (panel_logs.contains(name)) {
        return fps.collect { d -> ["logs/panel_resources/${name}${d.name}", d] }
    } else {
        return fps.collect { d -> ["logs/${meta.key}/${name}.${meta.id}${d.name}", d] }
    }

}
