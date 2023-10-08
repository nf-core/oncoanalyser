//
// Bam Tools calculates summary statistics for BAMs
//

import Constants
import Utils

include { BAMTOOLS } from '../../modules/local/bamtools/main'

workflow BAMTOOLS_METRICS {
    take:
        // Sample data
        ch_inputs      // channel: [mandatory] [ meta ]

        // Reference data
        genome_fasta   // channel: [mandatory] /path/to/genome_fasta
        genome_version // channel: [mandatory] genome version

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()


        // // NOTE(SW): I previously collapsed multiple normals to avoid processing duplication, however I'm now moving
        // // towards allowing multiple tumors per subject/group. This means that each normal should only be specified
        // // exactly once. I will not plan to support multiple tumor and normal combinations, this instead will required
        // // duplicate processing for some tasks but this workflow modality is likely to be a rare use-case.


        // Select inputs and then runnable/skip

        // Collapse by filepath
        //   * this must be generalised

        // Run process

        // Replicate by files
        //   * this must be generalised

        // Set all outputs i.e. add skip



        // Sort inputs
        // channel: runnable: [ meta ]
        // channel: skip: [ meta ]
        ch_inputs_sorted = ch_inputs
            .branch { meta ->

                def has_tumor_dna = Utils.hasTumorDnaBam(meta)
                def has_normal_dna = Utils.hasNormalDnaBam(meta)

                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAMTOOLS_TUMOR)

                runnable: (has_tumor_dna || has_normal_dna) && !has_existing
                skip: true
            }

        // Flatten into BAM/BAI pairs, select inputs that are eligible to run
        // channel: runnable: [ meta_extra, bam, bai ]
        // channel: skip: [ meta_extra ]
        ch_bams_bais_sorted = ch_inputs_sorted.runnable
            .flatMap { meta ->

                def tumor_bam = []
                def tumor_bai = []

                def normal_bam = []
                def normal_bai = []

                if (Utils.hasTumorDnaBam(meta)) {
                    tumor_bam = Utils.getTumorDnaBam(meta)
                    tumor_bai = Utils.getTumorDnaBai(meta)
                }

                if (Utils.hasNormalDnaBam(meta)) {
                    normal_bam = Utils.getNormalDnaBam(meta)
                    normal_bai = Utils.getNormalDnaBai(meta)
                }

                return [
                    [[key: meta.group_id, *:meta, sample_type: 'tumor'], tumor_bam, tumor_bai],
                    [[key: meta.group_id, *:meta, sample_type: 'normal'], normal_bam, normal_bai],
                ]
            }
            .branch { meta_extra, bam, bai ->
                runnable: bam && bai
                skip: true
                    return meta_extra
            }


        // Collapse duplicate files
        ch_bams_bais_dedup = ch_bams_bais_sorted.runnable
            .map { meta_extra, bam, bai ->
                return [[bam, bai], meta_extra]
            }
            .groupTuple()
            .map { fps, meta_extras ->

                def key = meta_extras.collect { meta_extra -> meta_extra.key }
                def meta_bamtools = [
                    id: key.join('__'),
                ]

                return [key, meta_extras, fps]
            }
            .multiMap { key, meta_extras, fps ->

            }
            .view()

        // Create process input channel






        /*




        // Select input sources
        // channel: [ meta_bamtools, bam, bai ]
        ch_bamtools_inputs_all = ch_inputs
            .flatMap { meta ->
                def inputs = []

                def meta_bamtools_tumor = [
                    key: meta.id,
                    id: Utils.getTumorDnaSampleName(meta),
                    // NOTE(SW): must use string representation for caching purposes
                    sample_type_str: Constants.SampleType.TUMOR.name(),
                ]
                def tumor_bam = Utils.getTumorDnaBam(meta)
                inputs.add([meta_bamtools_tumor, tumor_bam, "${tumor_bam}.bai"])

                if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
                    def meta_bamtools_normal = [
                        key: meta.id,
                        id: Utils.getNormalDnaSampleName(meta),
                        // NOTE(SW): must use string representation for caching purposes
                        sample_type_str: Constants.SampleType.NORMAL.name(),
                    ]
                    def normal_bam = Utils.getNormalDnaBam(meta)
                    inputs.add([meta_bamtools_normal, normal_bam, "${normal_bam}.bai"])
                }

                return inputs
            }

        // Collapse duplicate files e.g. repeated normal BAMs for multiple tumor samples
        // NOTE(SW): no effective blocking by .groupTuple() as we're not dependent
        // on any process
        // channel: [ meta_bamtools, bam, bai ]
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
            genome_fasta,
            genome_version,
        )

        // Set version
        ch_versions = ch_versions.mix(BAMTOOLS.out.versions)

        // Replicate outputs to reverse unique operation
        // channel: [ meta_bamtools_individual, sample_type_str, metrics ]
        ch_bamtools_out_individual = BAMTOOLS.out.metrics
            .flatMap { meta_bamtools_shared, metrics ->
                meta_bamtools_shared.keys.collect { key ->
                    return [meta_bamtools_shared + [key: key], meta_bamtools_shared.sample_type_str, metrics]
                }
            }

        // Set outputs
        ch_outputs = WorkflowOncoanalyser.restoreMeta(ch_bamtools_out_individual, ch_inputs)
            .branch { meta, sample_type_str, metrics ->
                def sample_type = Utils.getEnumFromString(sample_type_str, Constants.SampleType)
                somatic: sample_type == Constants.SampleType.TUMOR
                    return [meta, metrics]
                germline: sample_type == Constants.SampleType.NORMAL
                    return [meta, metrics]
            }



        */

        ch_outputs = ch_inputs
            .multiMap {
                somatic: it
                germline: it
            }




    emit:
        somatic  = ch_outputs.somatic  // channel: [ meta, metrics ]
        germline = ch_outputs.germline // channel: [ meta, metrics ]

        versions = ch_versions         // channel: [ versions.yml ]
}
