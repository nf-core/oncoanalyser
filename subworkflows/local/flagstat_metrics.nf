//
// SAMtools flagstat generates statistics for read alignments from the SAM FLAG field
//

import Constants
import Utils

include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'

workflow FLAGSTAT_METRICS {
    take:
        // Sample data
        ch_inputs  // channel: [mandatory] [ meta ]

        // Params
        run_config // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        /*

        // Select input source
        // Select input sources
        // channel: [ meta_flagstat, bam, bai ]
        ch_flagstat_inputs_all = ch_inputs
            .flatMap { meta ->
                def inputs = []

                def meta_flagstat_tumor = [
                    key: meta.id,
                    id: Utils.getTumorDnaSampleName(meta),
                    // NOTE(SW): must use string representation for caching purposes
                    sample_type_str: Constants.SampleType.TUMOR.name(),
                ]
                def tumor_bam = Utils.getTumorDnaBam(meta)
                inputs.add([meta_flagstat_tumor, tumor_bam, "${tumor_bam}.bai"])

                if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
                    def meta_flagstat_normal = [
                        key: meta.id,
                        id: Utils.getNormalDnaSampleName(meta),
                        // NOTE(SW): must use string representation for caching purposes
                        sample_type_str: Constants.SampleType.NORMAL.name(),
                    ]
                    def normal_bam = Utils.getNormalDnaBam(meta)
                    inputs.add([meta_flagstat_normal, normal_bam, "${normal_bam}.bai"])
                }

                return inputs
            }

        // Collapse duplicate files e.g. repeated normal BAMs for multiple tumor samples
        // channel: [ meta_flagstat_shared, bam, bai ]
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
        // channel: [ meta_flagstat_individual, flagstat ]
        ch_flagstat_out_individual = SAMTOOLS_FLAGSTAT.out.flagstat
            .flatMap { meta_flagstat_shared, flagstat ->
                def sample_type = Utils.getEnumFromString(meta_flagstat_shared.sample_type_str, Constants.SampleType)
                meta_flagstat_shared.keys.collect { key ->
                    return [meta_flagstat_shared + [key: key], sample_type, flagstat]
                }
            }

        // Set outputs
        // channel (somatic): [ meta, flagstat ]
        // channel (germline): [ meta, flagstat ]
        ch_outputs = WorkflowOncoanalyser.restoreMeta(ch_flagstat_out_individual, ch_inputs)
            .branch { meta, sample_type, flagstat ->
                somatic: sample_type == Constants.SampleType.TUMOR
                    return [meta, flagstat]
                germline: sample_type == Constants.SampleType.NORMAL
                    return [meta, flagstat]
            }
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

        */

        ch_outputs = ch_inputs
          .multiMap{ meta ->
              somatic: [meta, []]
              germline: [meta, []]
          }

    emit:
        somatic  = ch_outputs.somatic  // channel: [ meta, metrics ]
        germline = ch_outputs.germline // channel: [ meta, metrics ]

        versions = ch_versions         // channel: [ versions.yml ]
}
