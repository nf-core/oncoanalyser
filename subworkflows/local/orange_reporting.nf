//
// XXX
//
import Constants
import Utils

include { ORANGE } from '../../modules/local/orange/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'

workflow ORANGE_REPORTING {
    take:
        // Sample data
        ch_inputs
        ch_inputs_wgs
        ch_bamtools_somatic
        ch_bamtools_germline
        ch_chord
        ch_lilac
        ch_sage_somatic_tumor_bqr
        ch_sage_somatic_normal_bqr
        ch_sage_germline_coverage
        ch_purple
        ch_linx_somatic_annotation
        ch_linx_somatic_plot
        ch_linx_germline_annotation
        ch_protect
        ch_peach
        ch_cuppa
        ch_cuppa_feature_plot
        ch_cuppa_summary_plot
        ch_virusinterpreter

        // Reference data
        ref_data_genome_version
        ref_data_disease_ontology
        ref_data_known_fusion_data
        ref_data_driver_gene_panel
        ref_data_cohort_mapping
        ref_data_cohort_percentiles

        // Parameters
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // SAMtools flagstat
        // Select input source
        // channel (present): [val(meta), sample_type, flagstat]
        // channel (absent): [val(meta)]
        ch_inputs_flagstat = ch_inputs_wgs
            .flatMap { meta -> [Constants.SampleType.TUMOR, Constants.SampleType.NORMAL].collect { [meta, it] } }
            .branch { meta, sample_type ->
                def key = [Constants.FileType.FLAGSTAT, sample_type, Constants.SequenceType.WGS]
                present: meta.containsKey(key)
                    return [meta, sample_type, meta.getAt(key)]
                absent: ! meta.containsKey(key)
            }

        // Create inputs and create process-specific meta
        // channel: [val(meta_flagstat), bam, bai]
        ch_flagstat_inputs_all = ch_inputs_flagstat.absent
            .map { meta, sample_type ->
                def bam = meta.getAt([Constants.FileType.BAM, sample_type, Constants.SequenceType.WGS])
                def sample_name = meta.getAt(['sample_name', sample_type, Constants.SequenceType.WGS])
                def meta_flagstat = [
                    key: meta.id,
                    id: sample_name,
                    // NOTE(SW): must use string representation for caching purposes
                    sample_type_str: sample_type.name(),
                ]
                return [meta_flagstat, bam, "${bam}.bai"]
            }

        // Collapse duplicate files e.g. repeated normal BAMs for multiple tumor samples
        // channel: [val(meta_flagstat_shared), bam, bai]
        ch_flagstat_inputs = ch_flagstat_inputs_all
            .map { [it[1..-1], it[0]] }
            .groupTuple()
            .map { filepaths, meta_flagstat ->
                def (keys, sample_names, sample_type_strs) = meta_flagstat
                    .collect {
                        [it.key, it.id, it.sample_type_str]
                    }
                    .transpose()

                def sample_type_strs_unique = sample_type_strs.unique(false)
                assert sample_type_strs_unique.size() == 1
                def sample_type_str = sample_type_strs_unique[0]

                def meta_flagstat_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sample_type_str: sample_type_str,
                ]
                return [meta_flagstat_new, *filepaths]
            }

        // Run process
        SAMTOOLS_FLAGSTAT(
            ch_flagstat_inputs,
        )

        // Set version
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

        // Replicate outputs to reverse unique operation
        // channel: [val(meta_flagstat_individual), flagstat]
        ch_flagstat_out = SAMTOOLS_FLAGSTAT.out.flagstat
            .flatMap { meta_flagstat_shared, flagstat ->
                def sample_type = Utils.getEnumFromString(meta_flagstat_shared.sample_type_str, Constants.SampleType)
                meta_flagstat_shared.keys.collect { key ->
                    return [meta_flagstat_shared + [key: key], sample_type, flagstat]
                }
            }

        // Combine input flagstat channels, restoring original meta where required, split by sample type
        // channel (somatic): [val(meta), flagstat]
        // channel (germline): [val(meta), flagstat]
        ch_orange_inputs_flagstat = Channel.empty()
            .concat(
                ch_inputs_flagstat.present,
                WorkflowOncoanalyser.restoreMeta(ch_flagstat_out, ch_inputs),
            )
            .branch { meta, sample_type, flagstat ->
                somatic: sample_type == Constants.SampleType.TUMOR
                    return [meta, flagstat]
                germline: sample_type == Constants.SampleType.NORMAL
                    return [meta, flagstat]
            }

        // TODO(SW): handle CUPPA WTS, WGS, WGS

        // Select input source
        // NOTE(SW): we could consider not allowing inputs from the samplesheet here since this nothing follows
        ch_orange_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run.bamtools ? ch_bamtools_somatic : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.INPUT_BAMTOOLS_TXT_TUMOR),
            run.bamtools ? ch_bamtools_germline : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.INPUT_BAMTOOLS_TXT_NORMAL),
            ch_orange_inputs_flagstat.somatic,
            ch_orange_inputs_flagstat.germline,
            run.chord ? ch_chord : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CHORD),
            run.lilac ? ch_lilac : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LILAC_DIR),
            run.sage ? ch_sage_somatic_tumor_bqr : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_BQR_TUMOR),
            run.sage ? ch_sage_somatic_normal_bqr : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_BQR_NORMAL),
            run.sage ? ch_sage_germline_coverage : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_COVERAGE),
            run.purple ? ch_purple : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR),
            run.linx ? ch_linx_somatic_annotation : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            run.linx ? ch_linx_somatic_plot : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_PLOT_DIR_TUMOR),
            run.linx ? ch_linx_germline_annotation : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_PLOT_DIR_NORMAL),
            run.protect ? ch_protect : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PROTECT_TSV),
            run.peach ? ch_peach : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PEACH_TSV),
            run.cuppa ? ch_cuppa : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUTS_CUPPA_CSV),
            run.cuppa ? ch_cuppa_feature_plot : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CUPPA_FEATURE_PLOT),
            run.cuppa ? ch_cuppa_summary_plot : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CUPPA_SUMMARY_PLOT),
            run.virusinterpreter ? ch_virusinterpreter : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.VIRUSINTERPRETER_TSV),
        )

        ch_orange_inputs = ch_orange_inputs_source
            .map {
                def meta = it[0]
                def meta_orange = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorWgsSampleName(meta),
                    normal_id: Utils.getNormalWgsSampleName(meta),
                ]
                return [meta_orange, *it[1..-1]]
            }

        // Run process
        ORANGE(
            ch_orange_inputs,
            ref_data_genome_version,
            ref_data_disease_ontology,
            ref_data_known_fusion_data,
            ref_data_driver_gene_panel,
            ref_data_cohort_mapping,
            ref_data_cohort_percentiles,
            "5.31 [oncoanalyser]",
        )

        // Set outputs
        ch_versions = ch_versions.mix(ORANGE.out.versions)

    emit:
        versions  = ch_versions // channel: [versions.yml]
}
