import nextflow.Nextflow

class SampleSheet {

    public static parseInput(input_fp_str, stub_run, run_mode, log) {

        if (!input_fp_str) {
            log.error "Missing required --input argument"
            Nextflow.exit(1)
        }

        // NOTE(SW): using NF .splitCsv channel operator, hence should be easily interchangable with NF syntax

        def input_fp = getFileObject(input_fp_str)
        def inputs = nextflow.splitter.SplitterEx.splitCsv(input_fp, [header: true])
            .groupBy { it['group_id'] }
            .collect { group_id, entries ->

                def meta = [group_id: group_id]
                def sample_keys = [] as Set

                entries.each { entry ->
                    constructSampleMetaFromEntry(entry, meta, sample_keys, group_id, log)
                }

                // Checks per sample
                sample_keys.each { sample_key ->

                    if(stub_run) {
                        return
                    }

                    def meta_sample = meta[sample_key]

                    setCramPaths(meta_sample)
                    checkRawReadDataExists(meta_sample, group_id, log)
                    checkReduxTsvsExist(meta_sample, log)
                    checkAndSetFileIndexes(meta_sample, log)
                }

                // Checks per group
                disallowDuplicateSampleIds(meta, sample_keys, log)
                disallowInvalidSampleCombinations(meta, run_mode, log)
                checkLongitudinalSampleInputs(meta, log)

                return meta
            }

        return inputs
    }

    private static void constructSampleMetaFromEntry(entry, meta, sample_keys, group_id, log) {

        // Add subject id if absent or check if current matches existing
        if (meta.containsKey('subject_id') && meta.subject_id != entry.subject_id) {
            log.error "got unexpected subject name for ${group_id} ${meta.subject_id}: ${entry.subject_id}"
            Nextflow.exit(1)
        } else {
            meta.subject_id = entry.subject_id
        }

        def sample_type_enum = Enums.getValidatedEnumFromString(entry.sample_type, Constants.SampleType, log)
        def sequence_type_enum = Enums.getValidatedEnumFromString(entry.sequence_type, Constants.SequenceType, log)
        def filetype_enum = Enums.getValidatedEnumFromString(entry.filetype, Constants.FileType, log)

        def sample_key = [sample_type_enum, sequence_type_enum]
        sample_keys.add(sample_key) // Record sample key to simplify iteration later on

        def meta_sample = meta.get(sample_key, [:])

        // Info data
        def info_data = [:]
        if (entry.containsKey('info')) {
            // Parse
            entry.info
                .tokenize(';')
                .each { e ->
                    def (k, v) = e.tokenize(':')
                    def info_field_enum = Enums.getValidatedEnumFromString(k, Constants.InfoField, log)

                    if (info_data.containsKey(info_field_enum)) {
                        log.error "got duplicate info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${info_field_enum}"
                        Nextflow.exit(1)
                    }

                    if (!v && info_field_enum !== Constants.InfoField.LONGITUDINAL_SAMPLE) {
                        log.error "got empty value for ${group_id} ${sample_type_enum}/${sequence_type_enum} ${info_field_enum}"
                        Nextflow.exit(1)
                    }

                    info_data[info_field_enum] = v
                }

            // Process
            if (info_data.containsKey(Constants.InfoField.CANCER_TYPE)) {
                meta[Constants.InfoField.CANCER_TYPE] = info_data[Constants.InfoField.CANCER_TYPE]
            }

        }

        if (info_data.containsKey(Constants.InfoField.LONGITUDINAL_SAMPLE)) {

            if (meta_sample.containsKey('longitudinal_sample_id') && meta_sample.longitudinal_sample_id != entry.sample_id) {
                log.error "got multiple longitudinal samples for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${entry.sample_id}"
                Nextflow.exit(1)
            }

            meta_sample.longitudinal_sample_id = entry.sample_id

        } else if (meta_sample.containsKey('sample_id') && meta_sample.sample_id != entry.sample_id) {

            log.error "got unexpected sample name for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${entry.sample_id}"
            Nextflow.exit(1)

        } else {

            meta_sample.sample_id = entry.sample_id

        }

        // Filetype uniqueness
        if (meta_sample.containsKey(filetype_enum) & filetype_enum != Constants.FileType.FASTQ) {
            log.error "got duplicate file for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${filetype_enum}"
            Nextflow.exit(1)
        }

        // Handle inputs appropriately
        if (filetype_enum === Constants.FileType.FASTQ) {

            if (!info_data.containsKey(Constants.InfoField.LIBRARY_ID)) {
                log.error "missing 'library_id' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                Nextflow.exit(1)
            }

            if (!info_data.containsKey(Constants.InfoField.LANE)) {
                log.error "missing 'lane' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                Nextflow.exit(1)
            }

            def fastq_entries = entry.filepath.tokenize(';')

            if (fastq_entries.size() != 2) {
                log.error "expected exactly 2 FASTQ files delimited by ';' (i.e. '<fwd>;<rev>') but found ${fastq_entries.size} " +
                        " files for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                Nextflow.exit(1)
            }

            def (fwd, rev) = fastq_entries
            def fastq_key = [info_data[Constants.InfoField.LIBRARY_ID], info_data[Constants.InfoField.LANE]]

            if (!meta_sample.containsKey(filetype_enum)) {
                meta_sample[filetype_enum] = [:]
            }

            if (meta_sample[filetype_enum].containsKey(fastq_key)) {
                log.error "got duplicate lane + library_id data for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${fastq_key}"
                Nextflow.exit(1)
            }

            meta_sample[filetype_enum][fastq_key] = ['fwd': getFileObject(fwd), 'rev': getFileObject(rev)]

        } else {

            meta_sample[filetype_enum] = getFileObject(entry.filepath)

        }
    }

    private static getFileObject(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    private static void setCramPaths(meta_sample) {

        // CRAMs are passed to hmftools as if they were BAMs, e.g. `-bam_file /path/to/tumor.cram`
        // We therefore set the BAM/BAI path to be the CRAM/CRAI path

        if (meta_sample.containsKey(Constants.FileType.CRAM_REDUX)) {
            meta_sample[Constants.FileType.BAM_REDUX] = meta_sample.remove(Constants.FileType.CRAM_REDUX)
        }

        if (meta_sample.containsKey(Constants.FileType.CRAM)) {
            meta_sample[Constants.FileType.BAM] = meta_sample.remove(Constants.FileType.CRAM)
        }

        // The BAI key is used to store the index for both regular/REDUX CRAMs/BAMs
        if (meta_sample.containsKey(Constants.FileType.CRAI)) {
            meta_sample[Constants.FileType.BAI] = meta_sample.remove(Constants.FileType.CRAI)
        }

    }

    private static void checkAndSetFileIndexes(meta_sample, log) {

        meta_sample.keySet().each { key ->

            // NOTE(SW): I was going to use two maps but was unable to get an enum map to compile

            def index_enum
            def index_extension

            if (key === Constants.FileType.BAM || key === Constants.FileType.BAM_REDUX) {
                index_enum = Constants.FileType.BAI
                index_extension = 'bai'
            } else if (key === Constants.FileType.CRAM || key === Constants.FileType.CRAM_REDUX) {
                index_enum = Constants.FileType.CRAI
                index_extension = 'crai'
            } else if (key === Constants.FileType.ESVEE_VCF) {
                index_enum = Constants.FileType.ESVEE_VCF_TBI
                index_extension = 'tbi'
            } else if (key === Constants.FileType.SAGE_VCF) {
                index_enum = Constants.FileType.SAGE_VCF_TBI
                index_extension = 'tbi'
            } else {
                // Key not a file type, or not a file type that needs an index
                return
            }

            def index_already_provided = meta_sample.containsKey(index_enum)
            if (index_already_provided) {
                return
            }

            def file_path = meta_sample[key].toUriString()
            def index_path = nextflow.Nextflow.file("${file_path}.${index_extension}")

            if (!index_path.exists()) {
                log.error "no index provided or found for: ${file_path}"
                Nextflow.exit(1)
            }

            meta_sample[index_enum] = index_path
        }
    }

    private static void checkRawReadDataExists(meta_sample, group_id, log) {

        def missing_any_raw_read_data =
            !meta_sample.containsKey(Constants.FileType.BAM) &&
            !meta_sample.containsKey(Constants.FileType.BAM_REDUX) &&
            !meta_sample.containsKey(Constants.FileType.CRAM) &&
            !meta_sample.containsKey(Constants.FileType.CRAM_REDUX) &&
            !meta_sample.containsKey(Constants.FileType.FASTQ)

        if (missing_any_raw_read_data) {

            log.error "no BAM/CRAM nor BAM_REDUX/CRAM_REDUX nor FASTQ files provided for ${group_id} ${meta_sample.sample_id}\n\n" +
                    "NB: At least one of these files is required as they are the basis to determine input sample type."
            Nextflow.exit(1)
        }
    }

    private static void checkReduxTsvsExist(meta_sample, log) {

        if (!meta_sample.containsKey(Constants.FileType.BAM_REDUX)) {
            return
        }

        def bam_path = meta_sample[Constants.FileType.BAM_REDUX]
        def bam_dir = bam_path.getParent().toUriString()

        // Get user specified TSV paths
        def jitter_tsv = meta_sample[Constants.FileType.REDUX_JITTER_TSV]
        def ms_tsv = meta_sample[Constants.FileType.REDUX_MS_TSV]

        // If TSV paths not provided, default to TSV paths in the same dir as the BAM
        def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])
        jitter_tsv = jitter_tsv ?: "${bam_dir}/${sample_id}.jitter_params.tsv"
        ms_tsv = ms_tsv ?: "${bam_dir}/${sample_id}.ms_table.tsv.gz"

        jitter_tsv = nextflow.Nextflow.file(jitter_tsv)
        ms_tsv = nextflow.Nextflow.file(ms_tsv)

        def missing_tsvs = [:]
        if (!jitter_tsv.exists()) {
            missing_tsvs[Constants.FileType.REDUX_JITTER_TSV] = jitter_tsv
        }
        if (!ms_tsv.exists()) {
            missing_tsvs[Constants.FileType.REDUX_MS_TSV] = ms_tsv
        }

        if (missing_tsvs.size() > 0) {

            def error_message = []

            error_message.add("When only specifying filetype ${Constants.FileType.BAM_REDUX} in the sample sheet, make sure the REDUX BAM and TSVs are in the same dir:")
            error_message.add("${bam_path.toUriString()} (${Constants.FileType.BAM_REDUX})")
            missing_tsvs.each { error_message.add("${it.value} (missing expected ${it.key})") }
            error_message.add("")
            error_message.add(
                "Alternatively, provide the TSV paths in the sample sheet using filetype values: " +
                    "${Constants.FileType.REDUX_JITTER_TSV}, " +
                    "${Constants.FileType.REDUX_MS_TSV}"
            )

            log.error error_message.join("\n")
            Nextflow.exit(1)
        }

        // Set parsed REDUX TSV paths in metadata object
        meta_sample[Constants.FileType.REDUX_JITTER_TSV] = jitter_tsv
        meta_sample[Constants.FileType.REDUX_MS_TSV] = ms_tsv
    }

    private static void disallowDuplicateSampleIds(meta, sample_keys, log) {

        // Enforce unique samples names within groups
        def sample_ids_duplicated = sample_keys
            .groupBy { meta.getOrDefault(it, [:]).getOrDefault('sample_id', null) }
            .findResults { k, v -> k !== null & v.size() > 1 ? [k, v] : null }

        if (sample_ids_duplicated) {
            def duplicate_message_strs = sample_ids_duplicated.collect { sample_id, keys ->
                def key_strs = keys.collect { sample_type, sequence_type -> "${sample_type}/${sequence_type}" }
                return "  * ${sample_id}: ${key_strs.join(", ")}"
            }
            log.error "duplicate sample names found for ${meta.group_id}:\n\n${duplicate_message_strs.join("\n")}"
            Nextflow.exit(1)
        }
    }

    private static void disallowInvalidSampleCombinations(meta, run_mode, log) {

        // Do not allow normal DNA only
        if (Utils.hasNormalDna(meta) && !Utils.hasTumorDna(meta)) {
            log.error "found only normal DNA input for ${meta.group_id} but germline only analysis is not supported"
            Nextflow.exit(1)
        }

        // Do not allow CRAM RNA input
        if (Utils.hasTumorRnaBam(meta) && Utils.getTumorRnaBam(meta).toString().endsWith('cram')) {
            log.error "found tumor RNA CRAM input for ${meta.group_id} but RNA CRAM input is not supported"
            Nextflow.exit(1)
        }

        // Do not allow donor sample without normal sample
        if (Utils.hasDonorDna(meta) && !Utils.hasNormalDna(meta)) {
            log.error "a donor sample but not normal sample was found for ${meta.group_id}\n\n" +
                "Analysis with a donor sample requires a normal sample."
            Nextflow.exit(1)
        }

        // Apply some required restrictions to targeted mode
        if (run_mode === Constants.RunMode.TARGETED) {

            // Do not allow donor DNA
            if (Utils.hasDonorDna(meta)) {
                log.error "targeted mode is not compatible with the donor DNA BAM/CRAM provided for ${meta.group_id}\n\n" +
                    "The targeted workflow supports only tumor and normal DNA BAM/CRAMs (and tumor RNA BAM/CRAMs for TSO500)"
                Nextflow.exit(1)
            }

            // Do not allow only tumor RNA
            if (Utils.hasTumorRna(meta) && !Utils.hasTumorDna(meta)) {
                log.error "targeted mode is not compatible with only tumor RNA provided for ${meta.group_id}\n\n" +
                    "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                    "the configured panel."
                Nextflow.exit(1)
            }
        }
    }

    private static void checkLongitudinalSampleInputs(meta, log) {

        // For purity estimation with WISP, require primary normal DNA BAM when an AMBER directory is provided
        def meta_tumor_dna = meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
        def longitudinal = meta_tumor_dna.containsKey('longitudinal_sample_id')
        def has_amber_dir = meta_tumor_dna.containsKey(Constants.FileType.AMBER_DIR)
        def has_normal_dna_bam = Utils.hasNormalDnaBam(meta) || Utils.hasNormalDnaReduxBam(meta)

        if (longitudinal && has_amber_dir && !has_normal_dna_bam) {
            log.error "AMBER input was provided without the required primary normal DNA BAM for ${meta.group_id}"
            Nextflow.exit(1)
        }

    }
}
