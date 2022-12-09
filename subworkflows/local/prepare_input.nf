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
                    def sample_type_enum = Utils.getEnumFromString(it.sample_type, Constants.DataType)
                    if (!sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.DataType).join('\n  - ')
                        log.error "\nERROR: recieved invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
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
                    if (input_types_seen.contains([sample_type_enum, filetype_enum])) {
                        log.error "\nERROR: got duplicate inputs for ${key}: ${sample_type_enum}/${filetype_enum}"
                        System.exit(1)
                    }
                    input_types_seen.push([sample_type_enum, filetype_enum])

                    // Check for relevant indices
                    if (!workflow.stubRun) {
                        def filetype_bai = [
                            Constants.FileType.BAM_WGS,
                            Constants.FileType.BAM_WTS,
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
                        } else {
                          return
                        }

                        def index_fp_str = "${it.filepath}.${index_ext}".toString()
                        def index_fp = file(index_fp_str)
                        if (! index_fp.exists()) {
                            log.error "\nERROR: No index found for ${it.filepath}"
                            System.exit(1)
                        }
                    }

                    // Sample name
                    def key_sample_name = ['sample_name', sample_type_enum]
                    if (meta.containsKey(key_sample_name) && meta[key_sample_name] != it.sample_name) {
                        log.error "\nERROR: got unexpected sample name for ${key}/${meta[key_sample_name]}: ${it.sample_name}"
                        System.exit(1)
                    } else {
                        meta[key_sample_name] = it.sample_name
                    }

                    // Filepath
                    def key_file = [filetype_enum, sample_type_enum]
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
