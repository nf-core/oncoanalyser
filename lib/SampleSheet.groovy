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

                entries.each { entry ->
                    createOrUpdateSampleMeta(entry, meta, log)
                }

                // Checks per group
                disallowDuplicateSampleIds(meta, log)
                disallowInvalidSampleCombinations(meta, run_mode, log)
                checkLongitudinalSampleInputs(meta, log)

                return meta
            }

        return inputs
    }

    private static void createOrUpdateSampleMeta(entry, meta, log) {

        def group_id = meta.group_id

        // Add subject id if absent or check if current matches existing
        if (meta.containsKey('subject_id') && meta.subject_id != entry.subject_id) {
            log.error "got unexpected subject name for ${group_id} ${meta.subject_id}: ${entry.subject_id}"
            Nextflow.exit(1)
        } else {
            meta.subject_id = entry.subject_id
        }

        def sample_type = Enums.getValidatedEnumFromString(entry.sample_type, Constants.SampleType, log)
        def sequence_type = Enums.getValidatedEnumFromString(entry.sequence_type, Constants.SequenceType, log)
        def file_type = Enums.getValidatedEnumFromString(entry.filetype, Constants.FileType, log)

        def log_sample_id = "${group_id} ${sample_type}/${sequence_type}"

        def sample_key = [sample_type, sequence_type]
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
                        log.error "got duplicate info field for ${log_sample_id}: ${info_field_enum}"
                        Nextflow.exit(1)
                    }

                    if (!v && info_field_enum !== Constants.InfoField.LONGITUDINAL_SAMPLE) {
                        log.error "got empty value for ${log_sample_id} ${info_field_enum}"
                        Nextflow.exit(1)
                    }

                    info_data[info_field_enum] = v
                }

            // Process
            if (info_data.containsKey(Constants.InfoField.CANCER_TYPE)) {
                meta[Constants.InfoField.CANCER_TYPE] = info_data[Constants.InfoField.CANCER_TYPE]
            }

        }

        if (meta_sample.containsKey(file_type) & file_type != Constants.FileType.FASTQ) {
            log.error "got duplicate file for ${log_sample_id}: ${file_type}"
            Nextflow.exit(1)
        }

        if (file_type === Constants.FileType.FASTQ) {
            setFastqPaths(meta_sample, entry, info_data, log_sample_id, log)
        } else {
            meta_sample[file_type] = getFileObject(entry.filepath)
        }

        setCramPaths(meta_sample)

        setSampleIds(meta_sample, entry, info_data, log_sample_id, log)

        checkRawReadDataExists(meta_sample, group_id, log)

        checkReduxTsvsExist(meta_sample, log)

        checkAndSetFileIndexes(meta_sample, log)
    }

    private static getFileObject(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    private static void setSampleIds(meta_sample, entry, info_data, log_sample_id, log) {

        if (info_data.containsKey(Constants.InfoField.LONGITUDINAL_SAMPLE)) {

            if (meta_sample.containsKey('longitudinal_sample_id') && meta_sample.longitudinal_sample_id != entry.sample_id) {
                log.error "got multiple longitudinal samples for ${log_sample_id}: ${entry.sample_id}"
                Nextflow.exit(1)
            }

            meta_sample.longitudinal_sample_id = entry.sample_id

        } else if (meta_sample.containsKey('sample_id') && meta_sample.sample_id != entry.sample_id) {

            log.error "got unexpected sample name for ${log_sample_id}: ${entry.sample_id}"
            Nextflow.exit(1)

        } else {

            meta_sample.sample_id = entry.sample_id

        }
    }

    private static void setFastqPaths(meta_sample, entry, info_data, log_sample_id, log) {

        if (!info_data.containsKey(Constants.InfoField.LIBRARY_ID)) {
            log.error "missing 'library_id' info field for ${log_sample_id}"
            Nextflow.exit(1)
        }

        if (!info_data.containsKey(Constants.InfoField.LANE)) {
            log.error "missing 'lane' info field for ${log_sample_id}"
            Nextflow.exit(1)
        }

        def fastq_entries = entry.filepath.tokenize(';')

        if (fastq_entries.size() != 2) {
            log.error "expected exactly 2 FASTQ files delimited by ';' (i.e. '<fwd>;<rev>') but found ${fastq_entries.size} " +
                " files for ${log_sample_id}"
            Nextflow.exit(1)
        }

        def (fwd, rev) = fastq_entries
        def fastq_key = [info_data[Constants.InfoField.LIBRARY_ID], info_data[Constants.InfoField.LANE]]

        if (!meta_sample.containsKey(Constants.FileType.FASTQ)) {
            meta_sample[Constants.FileType.FASTQ] = [:]
        }

        if (meta_sample[Constants.FileType.FASTQ].containsKey(fastq_key)) {
            log.error "got duplicate lane + library_id data for ${log_sample_id}: ${fastq_key}"
            Nextflow.exit(1)
        }

        meta_sample[Constants.FileType.FASTQ][fastq_key] = ['fwd': getFileObject(fwd), 'rev': getFileObject(rev)]

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

        // NOTE(LN): Cast keys to list to avoid ConcurrentModificationException
        meta_sample.keySet().toList().each { key ->

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

        def jitter_tsv = meta_sample[Constants.FileType.REDUX_JITTER_TSV]
        def ms_tsv = meta_sample[Constants.FileType.REDUX_MS_TSV]
        def bqr_tsv = meta_sample[Constants.FileType.REDUX_BQR_TSV]

        // If TSV paths not provided, default to TSV paths in the same dir as the BAM
        def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

        jitter_tsv = jitter_tsv ?: nextflow.Nextflow.file("${bam_dir}/${sample_id}.redux.jitter_params.tsv")
        ms_tsv = ms_tsv ?: nextflow.Nextflow.file("${bam_dir}/${sample_id}.redux.ms_table.tsv.gz")
        bqr_tsv = bqr_tsv ?: nextflow.Nextflow.file("${bam_dir}/${sample_id}.redux.bqr.tsv")

        // Check for missing TSV files
        def missing_tsvs = [:]
        if (!bqr_tsv.exists()) {
            missing_tsvs[Constants.FileType.REDUX_BQR_TSV] = bqr_tsv
        }
        if (!jitter_tsv.exists()) {
            missing_tsvs[Constants.FileType.REDUX_JITTER_TSV] = jitter_tsv
        }
        if (!ms_tsv.exists()) {
            missing_tsvs[Constants.FileType.REDUX_MS_TSV] = ms_tsv
        }

        if (missing_tsvs.size() > 0) {

            def error_message = []

            missing_tsvs.each { error_message.add("Missing ${it.key}: ${it.value}") }
            error_message.add("")
            error_message.add("Make sure the REDUX BAM and TSVs are in the same dir, or explicitly provide")
            error_message.add("the TSV paths in the sample sheet using filetype values: ${missing_tsvs.keySet().join(", ")}")

            log.error error_message.join("\n")
            Nextflow.exit(1)
        }

        // Set parsed REDUX TSV paths in metadata object
        meta_sample[Constants.FileType.REDUX_BQR_TSV] = bqr_tsv
        meta_sample[Constants.FileType.REDUX_JITTER_TSV] = jitter_tsv
        meta_sample[Constants.FileType.REDUX_MS_TSV] = ms_tsv
    }

    private static void disallowDuplicateSampleIds(meta, log) {

        def sample_ids_seen = [] as Set
        def sample_ids_duplicated = [] as Set

        meta.each { key, maybe_meta_sample ->

            def is_meta_sample = (maybe_meta_sample instanceof Map) && maybe_meta_sample.containsKey('sample_id')

            if(!is_meta_sample) {
                return
            }

            def sample_id = maybe_meta_sample.sample_id

            if(sample_ids_seen.contains(sample_id)) {
                sample_ids_duplicated.add(sample_id)
            }

            sample_ids_seen.add(sample_id)
        }

        if (sample_ids_duplicated) {
            log.error "duplicate sample id(s) found for group_id(${meta.group_id}): ${sample_ids_duplicated.join(', ')}"
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
