//
// Prepare input data
//

import Constants
import Utils


workflow PREPARE_INPUT {
    take:
        ch_samplesheet // channel: [mandatory] /path/to/samplesheet

    main:
        ch_inputs = Channel.fromPath(ch_samplesheet)
            .splitCsv(header: true)
            .map { [it.group_id, it] }
            .groupTuple()
            .map { group_id, entries ->

                def meta = [group_id: group_id]
                def sample_keys = [] as Set

                // Process each entry
                entries.each {
                    // Add subject id if absent or check if current matches existing
                    if (meta.containsKey('subject_id') && meta.subject_id != it.subject_id) {
                        log.error "\nERROR: got unexpected subject name for ${group_id}/${meta.subject_id}: ${it.subject_id}"
                        System.exit(1)
                    } else {
                        meta.subject_id = it.subject_id
                    }


                    // Sample type
                    def sample_type_enum = Utils.getEnumFromString(it.sample_type, Constants.SampleType)
                    if (!sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.SampleType).join('\n  - ')
                        log.error "\nERROR: received invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        System.exit(1)
                    }

                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (!sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "\nERROR: received invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        System.exit(1)
                    }

                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (!filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "\nERROR: received invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        System.exit(1)
                    }

                    def sample_key = [sample_type_enum, sequence_type_enum]
                    def meta_sample = meta.get(sample_key, [sample_id: it.sample_id])

                    if (meta_sample.sample_id != it.sample_id) {
                        log.error "\nERROR: got unexpected sample name for ${group_id}/${sample_type_enum} (${sequence_type_enum}): ${sample_name}"
                        System.exit(1)
                    }

                    if (meta_sample.containsKey(filetype_enum)) {
                        log.error "\nERROR: got duplicate file for ${group_id}/${sample_type_enum} (${sequence_type_enum}): ${filetype_enum}"
                        System.exit(1)
                    }

                    meta_sample[filetype_enum] = it.filepath

                    // Record sample key to simplify iteration later on
                    sample_keys << sample_key

                }

                // Check that required indexes are provided or are accessible
                sample_keys.each { sample_key ->

                    if (workflow.stubRun) {
                        return
                    }

                    meta[sample_key]*.key.each { key ->

                        // NOTE(SW): I was going to use two maps but was unable to get an enum map to compile

                        def index_enum
                        def index_str

                        if (key === Constants.FileType.BAM) {
                            index_enum = Constants.FileType.BAI
                            index_str = 'bai'
                        } else if (key === Constants.FileType.GRIDSS_VCF) {
                            index_enum = Constants.FileType.GRIDSS_VCF_TBI
                            index_str = 'vcf'
                        } else if (key === Constants.FileType.GRIPSS_VCF) {
                            index_enum = Constants.FileType.GRIPSS_VCF_TBI
                            index_str = 'vcf'
                        } else if (key === Constants.FileType.GRIPSS_UNFILTERED_VCF) {
                            index_enum = Constants.FileType.GRIPSS_UNFILTERED_VCF_TBI
                            index_str = 'vcf'
                        } else {
                            return
                        }

                        if (meta[sample_key].containsKey(index_enum)) {
                            return
                        }

                        def fp = meta[sample_key][key]
                        def index_fp = file("${fp}.${index_str}")
                        if (!index_fp.exists()) {
                            log.error "\nERROR: No index provided or found for ${fp}"
                            System.exit(1)
                        }
                    }
                }

                return meta
            }

    emit:
      data = ch_inputs
}
