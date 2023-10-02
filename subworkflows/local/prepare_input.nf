//
// Prepare input data
//

import Constants
import Utils


workflow PREPARE_INPUT {
    take:
        // Sample data
        ch_samplesheet // channel: [mandatory] /path/to/samplesheet

        // Params
        run_config     // channel: [mandatory] run configuration

    main:

        if (run_config.type == Constants.RunType.TUMOR_ONLY) {
            sample_types_allowed = [Constants.SampleType.TUMOR]
        } else if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
            sample_types_allowed = [
                Constants.SampleType.TUMOR,
                Constants.SampleType.NORMAL,
                Constants.SampleType.TUMOR_NORMAL,
            ]
        } else {
            assert false
        }

        if (run_config.mode == Constants.RunMode.DNA) {
            sequence_types_allowed = [Constants.SequenceType.DNA]
        } else if (run_config.mode == Constants.RunMode.RNA) {
            sequence_types_allowed = [Constants.SequenceType.RNA]
        } else if (run_config.mode == Constants.RunMode.DNA_RNA) {
            sequence_types_allowed = [Constants.SequenceType.DNA, Constants.SequenceType.RNA]
        } else {
            assert false
        }




        // TODO(SW): ensure that config enforces targeted cannot be RNA only




        ch_inputs = Channel.of(ch_samplesheet)
            .splitCsv(header: true)
            .map { [it.id, it] }
            .groupTuple()
            .map { key, entries ->
                def meta = [id: key]
                def input_types_seen = []

                // Process each entry
                entries.each {
                    // Add subject name if absent or check if current matches existing
                    if (meta.containsKey('subject_name') && meta.subject_name != it.subject_name) {
                        log.error "\nERROR: got unexpected subject name for ${key}/${meta.subject_name}: ${it.subject_name}"
                        System.exit(1)
                    } else {
                        meta.subject_name = it.subject_name
                    }


                    // Sample type
                    def sample_type_enum = Utils.getEnumFromString(it.sample_type, Constants.SampleType)
                    if (!sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.SampleType).join('\n  - ')
                        log.error "\nERROR: received invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        System.exit(1)
                    }

                    if (!sample_types_allowed.contains(sample_type_enum)) {
                        def sample_type_str = sample_types_allowed.collect { it.name().toLowerCase() }.join('\n  - ')
                        def run_type_str = run_config.type.name().toLowerCase()
                        log.error "\nERROR: received invalid sample type for ${run_type_str} " +
                            "run: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        System.exit(1)
                    }




                    // TODO(SW): understand whether this does (or why it doesn't) failure on DNA_RNA inputs (e.g. CUPPA)




                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (!sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "\nERROR: received invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        System.exit(1)
                    }

                    if (!sequence_types_allowed.contains(sequence_type_enum)) {
                        def sequence_types_str = sequence_types_allowed.collect { it.name().toLowerCase() }.join('\n  - ')
                        def run_mode_str = run_config.mode.name().toLowerCase()
                        log.error "\nERROR: received invalid sample mode for ${run_mode_str} " +
                            "run: '${it.sequence_type}'. Valid options are:\n  - ${sequence_types_str}"
                        System.exit(1)
                    }


                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (!filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "\nERROR: received invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        System.exit(1)
                    }

                    // Check whether this input type already exists
                    def key_input_type = [sample_type_enum, sequence_type_enum, filetype_enum]
                    if (input_types_seen.contains(key_input_type)) {
                        log.error "\nERROR: got duplicate inputs for ${key}: ${sample_type_enum}/${filetype_enum}"
                        System.exit(1)
                    }
                    input_types_seen.push(key_input_type)




                    // TODO(SW): implement index finding when not provided in samplesheet




                    // Check for relevant indices
                    if (!workflow.stubRun) {
                        def filetype_bai = [
                            Constants.FileType.BAM,
                        ]

                        def filetype_tbi = [
                            Constants.FileType.GRIDSS_VCF,
                            Constants.FileType.GRIPSS_VCF,
                            Constants.FileType.GRIPSS_UNFILTERED_VCF,
                        ]

                        def index_ext
                        if (filetype_bai.contains(filetype_enum)) {
                            index_ext = 'bai'
                        } else if (filetype_tbi.contains(filetype_enum)) {
                            index_ext = 'tbi'
                        }

                        if (index_ext) {
                            def index_fp_str = "${it.filepath}.${index_ext}".toString()
                            def index_fp = file(index_fp_str)
                            if (! index_fp.exists()) {
                                log.error "\nERROR: No index found for ${it.filepath}"
                                System.exit(1)
                            }
                        }
                    }

                    // Sample name
                    def key_sample_name
                    if (sample_type_enum == Constants.SampleType.TUMOR_NORMAL) {
                        if (it.sample_name.contains(';')) {
                            def sample_name_tokens = it.sample_name.split(';')
                            if (sample_name_tokens.size() != 2) {
                                log.error "\nERROR: expected two sample names but got ${sample_name_tokens.size()}"
                                System.exit(1)
                            }

                            def sample_name_and_types = [
                                [Constants.SampleType.TUMOR, Constants.SampleType.NORMAL],
                                sample_name_tokens,
                            ]
                            sample_name_and_types
                                .transpose()
                                .each { st, sn ->
                                    key_sample_name = ['sample_name', st, sequence_type_enum]
                                    meta[key_sample_name] = process_sample_name(sn, key_sample_name, meta)
                                }
                        } else {
                            key_sample_name = ['sample_name', Constants.SampleType.TUMOR, sequence_type_enum]
                            meta[key_sample_name] = process_sample_name(it.sample_name, key_sample_name, meta)
                        }
                    } else {
                        key_sample_name = ['sample_name', sample_type_enum, sequence_type_enum]
                        meta[key_sample_name] = process_sample_name(it.sample_name, key_sample_name, meta)
                    }


                    // Filepath
                    def key_file = [filetype_enum, sample_type_enum, sequence_type_enum]
                    if (meta.containsKey(key_file)) {
                        log.error "\nERROR: got duplicate file for ${key}: ${filetype_enum}/${sample_type_enum}"
                        System.exit(1)
                    } else {
                        meta[key_file] = it.filepath
                    }
                }
                return meta
            }




            /*
            // Check we have the required sample types for the specified run configuration
            ch_inputs
                .map { meta ->

                    def sample_types = []
                    for (e in meta) {
                        if (! e.key instanceof List || e.key[0] != 'sample_name') {
                            continue
                        }
                        sample_types.add(e.key[1..-1])
                    }

                    def required_sample_types
                    if (run_config.mode == Constants.RunMode.PANEL) {

                        required_sample_types = [
                            [Constants.SampleType.TUMOR, Constants.SequenceType.TARGETTED],
                        ]

                    } else if (run_config.mode == Constants.RunMode.DNA) {

                        if (run_config.type == Constants.RunType.TUMOR_ONLY) {
                            required_sample_types = [
                                [Constants.SampleType.TUMOR, Constants.SequenceType.DNA],
                            ]
                        } else if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
                            required_sample_types = [
                                [Constants.SampleType.TUMOR, Constants.SequenceType.DNA],
                                [Constants.SampleType.NORMAL, Constants.SequenceType.DNA],
                            ]
                        } else {
                            assert false
                        }

                    } else if (run_config.mode == Constants.RunMode.RNA) {

                        required_sample_types = [
                            [Constants.SampleType.TUMOR, Constants.SequenceType.RNA],
                        ]

                    } else if (run_config.mode == Constants.RunMode.DNA_RNA) {

                        if (run_config.type == Constants.RunType.TUMOR_ONLY) {
                            required_sample_types = [
                                [Constants.SampleType.TUMOR, Constants.SequenceType.DNA],
                                [Constants.SampleType.TUMOR, Constants.SequenceType.RNA],
                            ]
                        } else if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
                            required_sample_types = [
                                [Constants.SampleType.TUMOR, Constants.SequenceType.DNA],
                                [Constants.SampleType.TUMOR, Constants.SequenceType.RNA],
                                [Constants.SampleType.NORMAL, Constants.SequenceType.DNA],
                            ]
                        } else {
                            assert false
                        }

                    } else {
                        assert false
                    }

                    def sample_types_missing = required_sample_types - sample_types
                    def sample_types_extra = sample_types - required_sample_types


                    if (sample_types_missing) {

                        def sample_type_str = sample_types_missing
                            .collect { e ->
                                def (sample, sequence) = e.collect { it.name().toLowerCase() }
                                return "${sample}/${sequence}"
                            }
                            .join('\n  - ')

                        def run_type_str = run_config.type.name().toLowerCase()
                        log.error "\nERROR: missing required input for ${run_type_str} run:\n  - ${sample_type_str}"
                        System.exit(1)
                    }

                    // NOTE(SW): this shold never evalutate as true with the above checks in place
                    if (sample_types_extra) {

                        def sample_type_str = sample_types_extra
                            .collect { e ->
                                def (sample, sequence) = e.collect { it.name().toLowerCase() }
                                return "${sample}/${sequence}"
                            }
                            .join('\n  - ')

                        def run_type_str = run_config.type.name().toLowerCase()
                        log.error "\nERROR: extra input for ${run_type_str} run found:\n  - ${sample_type_str}"
                        System.exit(1)
                    }

                }
            */




    emit:
      data = ch_inputs
}

def process_sample_name(sample_name, key_sample_name, meta) {
    if (meta.containsKey(key_sample_name) && meta[key_sample_name] != sample_name) {
        log.error "\nERROR: got unexpected sample name for ${key}/${meta[key_sample_name]}: ${sample_name}"
        System.exit(1)
    }
    return sample_name
}
