//
// XXX
//
import Constants
import Utils

include { BAMTOOLS } from '../../modules/local/bamtools/main'

workflow BAMTOOLS_METRICS {
    take:
        // Sample data
        ch_inputs

        // Reference data
        ref_data_genome_fasta
        ref_data_genome_version

        // Params
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // NOTE(SW): CUPPA only requires metrics for the tumor sample in the upstream
        // process Virus Interpreter but ORANGE currently requires metrics for both tumor
        // and normal sample
        // channel: [val(meta_bamtools), bam, bai]
        ch_bamtools_inputs_all = ch_inputs
            .flatMap { meta ->
                def sample_types
                if (run.orange) {
                    sample_types = [Constants.SampleType.TUMOR, Constants.SampleType.NORMAL]
                } else {
                    sample_types = [Constants.SampleType.TUMOR]
                }

                return sample_types
                    .collect { sample_type ->
                        def bam = meta.get([Constants.FileType.BAM, sample_type, Constants.SequenceType.WGS])
                        def sample_name = meta.get(['sample_name', sample_type, Constants.SequenceType.WGS])
                        def meta_bamtools = [
                            key: meta.id,
                            id: sample_name,
                            // NOTE(SW): must use string representation for caching purposes
                            sample_type_str: sample_type.name(),
                        ]
                        return [meta_bamtools, bam, "${bam}.bai"]
                    }
            }

        // Collapse duplicate files e.g. repeated normal BAMs for multiple tumor samples
        // NOTE(SW): no effective blocking by .groupTuple() as we're not dependent
        // on any process
        // channel: [val(meta_bamtools), bam, bai]
        ch_bamtools_inputs = ch_bamtools_inputs_all
            .map { [it[1..-1], it[0]] }
            .groupTuple()
            .map { filepaths, meta_bamtools ->
                def (keys, sample_names, sample_type_strs) = meta_bamtools
                    .collect {
                        [it.key, it.id, it.sample_type_str]
                    }
                    .transpose()

                def sample_type_strs_unique = sample_type_strs.unique(false)
                assert sample_type_strs_unique.size() == 1
                def sample_type_str = sample_type_strs_unique[0]

                def meta_bamtools_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sample_type_str: sample_type_str,
                ]
                return [meta_bamtools_new, *filepaths]
            }

        // Run process
        BAMTOOLS(
            ch_bamtools_inputs,
            ref_data_genome_fasta,
            ref_data_genome_version,
        )

        // Set outputs, process outputs and restore original meta
        ch_versions = ch_versions.mix(BAMTOOLS.out.versions)

        // Replicate outputs to reverse unique operation
        // channel: [val(meta_bamtools_individual), sample_type_str, metrics]
        ch_bamtools_out_individual = BAMTOOLS.out.metrics
            .flatMap { meta_bamtools_shared, metrics ->
                meta_bamtools_shared.keys.collect { key ->
                    return [meta_bamtools_shared + [key: key], meta_bamtools_shared.sample_type_str, metrics]
                }
            }

        // Set outputs
        // channel (somatic): [val(meta), metrics]
        // channel (germline): [val(meta), metrics]
        ch_outputs = WorkflowOncoanalyser.restoreMeta(ch_bamtools_out_individual, ch_inputs)
            .branch { meta, sample_type_str, metrics ->
                def sample_type = Utils.getEnumFromString(sample_type_str, Constants.SampleType)
                somatic: sample_type == Constants.SampleType.TUMOR
                    return [meta, metrics]
                germline: sample_type == Constants.SampleType.NORMAL
                    return [meta, metrics]
            }

    emit:
        somatic  = ch_outputs.somatic  // channel: [val(meta), metrics]
        germline = ch_outputs.germline // channel: [val(meta), metrics]

        versions = ch_versions         // channel: [versions.yml]
}
