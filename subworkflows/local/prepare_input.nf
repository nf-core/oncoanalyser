import Constants
import Utils


workflow PREPARE_INPUT {
    take:
        ch_samplesheet

    main:
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
                        log.error "\nERROR: recieved invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        System.exit(1)
                    }

                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (!sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "\nERROR: recieved invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        System.exit(1)
                    }

                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (!filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "\nERROR: recieved invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        System.exit(1)
                    }

                    // Check whether this input type already exists
                    def key_input_type = [sample_type_enum, sequence_type_enum, filetype_enum]
                    if (input_types_seen.contains(key_input_type)) {
                        log.error "\nERROR: got duplicate inputs for ${key}: ${sample_type_enum}/${filetype_enum}"
                        System.exit(1)
                    }
                    input_types_seen.push(key_input_type)

                    // Check for relevant indices
                    if (!workflow.stubRun) {
                        def filetype_bai = [
                            Constants.FileType.BAM,
                        ]

                        def filetype_tbi = [
                            Constants.FileType.GRIDSS_VCF,
                            Constants.FileType.GRIPSS_HARD_VCF,
                            Constants.FileType.GRIPSS_SOFT_VCF,
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
