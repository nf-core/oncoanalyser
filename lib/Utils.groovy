//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

import nextflow.Nextflow
import nextflow.splitter.SplitterEx

class Utils {

    public static parseInput(input_fp_str, stub_run, log) {

        // NOTE(SW): using Nextflow .splitCsv channel operator, hence sould be easily interchangable

        def input_fp = nextflow.Nextflow.file(input_fp_str)
        def inputs = nextflow.splitter.SplitterEx.splitCsv(input_fp, [header: true])
            .groupBy { it['group_id'] }
            .collect { group_id, entries ->

                def meta = [group_id: group_id]
                def sample_keys = [] as Set

                // Process each entry
                entries.each {
                    // Add subject id if absent or check if current matches existing
                    if (meta.containsKey('subject_id') && meta.subject_id != it.subject_id) {
                        log.error "got unexpected subject name for ${group_id} ${meta.subject_id}: ${it.subject_id}"
                        System.exit(1)
                    } else {
                        meta.subject_id = it.subject_id
                    }


                    // Sample type
                    def sample_type_enum = Utils.getEnumFromString(it.sample_type, Constants.SampleType)
                    if (!sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.SampleType).join('\n  - ')
                        log.error "received invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        System.exit(1)
                    }

                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (!sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "received invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        System.exit(1)
                    }

                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (!filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "received invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        System.exit(1)
                    }

                    def sample_key = [sample_type_enum, sequence_type_enum]
                    def meta_sample = meta.get(sample_key, [sample_id: it.sample_id])

                    if (meta_sample.sample_id != it.sample_id) {
                        log.error "got unexpected sample name for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${it.sample_id}"
                        System.exit(1)
                    }

                    if (meta_sample.containsKey(filetype_enum)) {
                        log.error "got duplicate file for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${filetype_enum}"
                        System.exit(1)
                    }

                    meta_sample[filetype_enum] = Utils.getFileObject(it.filepath)

                    // Record sample key to simplify iteration later on
                    sample_keys << sample_key
                }

                // Check that required indexes are provided or are accessible
                sample_keys.each { sample_key ->

                    meta[sample_key]*.key.each { key ->

                        // NOTE(SW): I was going to use two maps but was unable to get an enum map to compile

                        def index_enum
                        def index_str

                        if (key === Constants.FileType.BAM) {
                            index_enum = Constants.FileType.BAI
                            index_str = 'bai'
                        } else if (key === Constants.FileType.GRIDSS_VCF) {
                            index_enum = Constants.FileType.GRIDSS_VCF_TBI
                            index_str = 'tbi'
                        } else if (key === Constants.FileType.GRIPSS_VCF) {
                            index_enum = Constants.FileType.GRIPSS_VCF_TBI
                            index_str = 'tbi'
                        } else if (key === Constants.FileType.GRIPSS_UNFILTERED_VCF) {
                            index_enum = Constants.FileType.GRIPSS_UNFILTERED_VCF_TBI
                            index_str = 'tbi'
                        } else if (key === Constants.FileType.SAGE_VCF) {
                            index_enum = Constants.FileType.SAGE_VCF_TBI
                            index_str = 'tbi'
                        } else {
                            return
                        }

                        if (meta[sample_key].containsKey(index_enum)) {
                            return
                        }

                        def fp = meta[sample_key][key].toUriString()
                        def index_fp = nextflow.Nextflow.file("${fp}.${index_str}")

                        if (!index_fp.exists() && !stub_run) {
                            def (sample_type, sequence_type) = sample_key
                            log.error "no index provided or found for ${meta.group_id} ${sample_type}/${sequence_type}: ${key}: ${fp}"
                            System.exit(1)
                        }

                        meta[sample_key][index_enum] = index_fp

                    }
                }

                return meta
            }

        return inputs
    }

    public static void validateInput(inputs, run_config, log) {

        def sample_keys = [
            [Constants.SampleType.TUMOR, Constants.SequenceType.DNA],
            [Constants.SampleType.TUMOR, Constants.SequenceType.RNA],
            [Constants.SampleType.NORMAL, Constants.SequenceType.DNA],
        ]

        inputs.each { meta ->

            // Require BAMs for each defined sample type
            // NOTE(SW): repeating key pairs above to avoid having to duplicate error messages
            sample_keys.each { key ->

                if (!meta.containsKey(key)) {
                    return
                }

                def (sample_type, sequence_type) = key

                if (!meta[key].containsKey(Constants.FileType.BAM)) {
                    log.error "no BAM provided for ${meta.group_id} ${sample_type}/${sequence_type}\n\n" +
                        "NB: BAMs are always required as they are the basis to determine input sample type."
                    System.exit(1)
                }

            }

            // Apply some required restrictions to targeted mode
            if (run_config.mode === Constants.RunMode.TARGETED) {

                // Do not allow normal DNA
                if (Utils.hasNormalDnaBam(meta)) {
                    log.error "targeted mode is not compatible with the normal DNA BAM provided for ${meta.group_id}\n\n" +
                        "The targeted workflow supports only tumor DNA BAMs (and tumor RNA BAMs for TSO500)"
                    System.exit(1)
                }

                // Do not allow only tumor RNA
                if (Utils.hasTumorRnaBam(meta) && !Utils.hasTumorDnaBam(meta)) {
                    log.error "targeted mode is not compatible with only tumor RNA provided for ${meta.group_id}\n\n" +
                        "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                        "the configured panel."
                    System.exit(1)
                }

                // Restrict tumor RNA inputs to the TSO500 panel
                if (Utils.hasTumorRnaBam(meta) && run_config.panel != 'tso500') {
                    def panel = run_config.panel.toUpperCase()
                        "Only the TSO500 panel supports tumor RNA analysis"
                    System.exit(1)
                }

            }

            // Do not allow normal DNA only
            if (Utils.hasNormalDnaBam(meta) && !Utils.hasTumorDnaBam(meta)) {
                log.error "germline only mode not supported, found only a normal DNA BAM for ${meta.group_id}\n"
                System.exit(1)
            }

            // Enforce unique samples names within groups
            def sample_ids_duplicated = sample_keys
                .groupBy { meta.getOrDefault(it, [:]).getOrDefault('sample_id', null) }
                .findResults { k, v -> k !== null & v.size() > 1 ? [k, v] : null }

            if (sample_ids_duplicated) {
                def duplicate_message_strs = sample_ids_duplicated.collect { sample_id, keys ->
                    def key_strs = keys.collect { sample_type, sequence_type -> "${sample_type}/${sequence_type}" }
                    return "  * ${sample_id}: ${key_strs.join(", ")}"
                }
                log.error "duplicate sample names found for ${meta.group_id}:\n\n" +
                    "${duplicate_message_strs.join("\n")}"
                System.exit(1)
            }

        }
    }

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    static public getEnumFromString(s, e) {
        try {
            return e.valueOf(s.toUpperCase())
        } catch(java.lang.IllegalArgumentException err) {
            return null
        }
    }

    public static getEnumNames(e) {
        e
            .values()
            *.name()
            *.toLowerCase()
    }


    static public getFileObject(path) {
        return path ? Nextflow.file(path) : []
    }


    // Sample names
    static public getTumorDnaSampleName(meta) {
        def meta_sample = meta[Constants.SampleType.TUMOR, Constants.SequenceType.DNA]
        return meta_sample['sample_id']
    }

    static public getTumorRnaSampleName(meta) {
        def meta_sample = meta[Constants.SampleType.TUMOR, Constants.SequenceType.RNA]
        return meta_sample['sample_id']
    }

    static public getNormalDnaSampleName(meta) {
        def meta_sample = meta[Constants.SampleType.NORMAL, Constants.SequenceType.DNA]
        return meta_sample['sample_id']
    }


    // Files
    static public getTumorDnaBam(meta) {
        def meta_sample = meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
        return meta_sample.getOrDefault(Constants.FileType.BAM, null)
    }

    static public getTumorDnaBai(meta) {
        def meta_sample = meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
        return meta_sample.getOrDefault(Constants.FileType.BAI, null)
    }

    static public hasTumorDnaBam(meta) {
        return getTumorDnaBam(meta) !== null
    }

    static public getTumorRnaBam(meta) {
        def meta_sample = meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.RNA], [:])
        return meta_sample.getOrDefault(Constants.FileType.BAM, null)
    }

    static public getTumorRnaBai(meta) {
        def meta_sample = meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.RNA], [:])
        return meta_sample.getOrDefault(Constants.FileType.BAI, null)
    }

    static public hasTumorRnaBam(meta) {
        return getTumorRnaBam(meta) !== null
    }


    static public getNormalDnaBam(meta) {
        def meta_sample = meta.getOrDefault([Constants.SampleType.NORMAL, Constants.SequenceType.DNA], [:])
        return meta_sample.getOrDefault(Constants.FileType.BAM, null)
    }

    static public getNormalDnaBai(meta) {
        def meta_sample = meta.getOrDefault([Constants.SampleType.NORMAL, Constants.SequenceType.DNA], [:])
        return meta_sample.getOrDefault(Constants.FileType.BAI, null)
    }

    static public hasNormalDnaBam(meta) {
        return getNormalDnaBam(meta) !== null
    }


    static public getRunMode(run_mode, log) {
        def run_mode_enum = Utils.getEnumFromString(run_mode, Constants.RunMode)
        if (!run_mode_enum) {
            def run_modes_str = Utils.getEnumNames(Constants.RunMode).join('\n  - ')
            log.error "recieved an invalid run mode: '${run_mode}'. Valid options are:\n  - ${run_modes_str}"
            System.exit(1)
        }
        return run_mode_enum
    }



    public static getInput(meta, key) {

        def result
        def (key_filetype, key_filetypes, key_sequencetypes) = key

        for (key_sample in [key_filetypes, key_sequencetypes].combinations()) {
            if (meta.containsKey(key_sample) && meta[key_sample].containsKey(key_filetype)) {
                return meta[key_sample].getAt(key_filetype)
            }
        }
    }

    public static hasExistingInput(meta, key) {
        return getInput(meta, key) !== null
    }

    public static selectCurrentOrExisting(val, meta, key) {
        if (hasExistingInput(meta, key)) {
          return getInput(meta, key)
        } else {
          return val
        }
    }

}
