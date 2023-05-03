//
// XXX
//
import Constants
import Utils

include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'

workflow FLAGSTAT_METRICS {
    take:
        ch_inputs

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input source
        // channel (present): [val(meta), sample_type, flagstat]
        // channel (absent): [val(meta)]
        ch_inputs_flagstat = ch_inputs
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
        ch_outputs = Channel.empty()
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
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    emit:
        somatic  = ch_outputs.somatic  // channel: [val(meta), metrics]
        germline = ch_outputs.germline // channel: [val(meta), metrics]

        versions = ch_versions         // channel: [versions.yml]
}
