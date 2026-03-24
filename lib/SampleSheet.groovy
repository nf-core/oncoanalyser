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

                /*
                // Example typical structure of `meta` once fully constructed
                // Contains the run metadata for every sample within the same group_id

                [
                    group_id:<string>,
                    subject_id:<string>,

                    // `meta_sample` (multiple entries shown)
                    [TUMOR, DNA]:[sample_id:<string>, BAM:<path_as_string>, BAI:<path_as_string>],
                    [NORMAL, DNA]:[sample_id:<string>, BAM:<path_as_string>, BAI:<path_as_string>],
                    [TUMOR, RNA]:[sample_id:<string>, FASTQ:[[S1, 001]:[fwd:<path_as_string>, rev:<path_as_string>]]]
                ]
                */

                entries.each { entry ->
                    createOrUpdateSampleMeta(entry, meta, sample_keys, log)
                }

                // Per sample checks once meta_sample objects are fully constructed
                sample_keys.each { sample_key ->

                    if(stub_run) {
                        // NOTE(LN): TODO: not sure if this skip is required
                        return
                    }

                    def meta_sample = meta[sample_key]

                    // NOTE(LN): The order in which these methods are executed is important
                    checkAndSetFileIndexes(meta_sample, log)
                    setCramPaths(meta_sample)
                    checkRawReadDataExists(meta_sample, group_id, log)
                    setReduxTsvDirIfUnset(meta_sample)

                }

                // Per group checks
                disallowDuplicateSampleIds(meta, log)
                disallowInvalidSampleCombinations(meta, run_mode, log)
                checkLongitudinalSampleInputs(meta, log)

                return meta
            }

        return inputs
    }

    private static void createOrUpdateSampleMeta(entry, meta, sample_keys, log) {

        def group_id = meta.group_id

        def sample_type = Enums.getValidatedEnumFromString(entry.sample_type, SampleMeta.SampleType, log)
        def sequence_type = Enums.getValidatedEnumFromString(entry.sequence_type, SampleMeta.SequenceType, log)
        def file_type = Enums.getValidatedEnumFromString(entry.filetype, SampleMeta.FileType, log)

        def log_sample_id = "group_id(${group_id}) sample_id(${sample_type})"
        def log_group_id = "group_id(${group_id})"

        setAndCheckSubjectId(meta, entry, log_group_id, log)

        def sample_key = [sample_type, sequence_type]
        sample_keys.add(sample_key) // Store to iterate over later

        def meta_sample = meta.get(sample_key, [:])

        // Handle info field
        def info_data = parseInfoField(entry, log_sample_id, log)

        if (info_data.containsKey(SampleMeta.InfoField.CANCER_TYPE)) {
            meta[SampleMeta.InfoField.CANCER_TYPE] = info_data[SampleMeta.InfoField.CANCER_TYPE]
        }

        setSampleIds(meta_sample, entry, info_data, log_group_id, log)

        // Set file paths
        if (meta_sample.containsKey(file_type) & file_type != SampleMeta.FileType.FASTQ) {
            log.error "got duplicate filetype(${file_type}) for ${log_sample_id}"
            Nextflow.exit(1)
        }

        if (file_type === SampleMeta.FileType.FASTQ) {
            setFastqPaths(meta_sample, entry, info_data, log_sample_id, log)
        } else {
            meta_sample[file_type] = getFileObject(entry.filepath)
        }
    }

    private static void setAndCheckSubjectId(meta, entry, log_group_id, log) {

        if(!meta.containsKey('subject_id')) {
            meta.subject_id = entry.subject_id
        }

        if (meta.subject_id != entry.subject_id) {
            log.error "got unexpected subject_id(${entry.subject_id}) for ${log_group_id}. " +
                "All samples within a group must have the same subject_id"

            Nextflow.exit(1)
        }
    }

    private static Map parseInfoField(entry, log_sample_id, log) {
        def info_data = [:]

        if (!entry.containsKey('info')) {
            return info_data
        }

        entry.info
            .tokenize(';')
            .each { e ->
                def (k, v) = e.tokenize(':')
                def info_field_enum = Enums.getValidatedEnumFromString(k, SampleMeta.InfoField, log)

                if (info_data.containsKey(info_field_enum)) {
                    log.error "got duplicate info field(${info_field_enum}) for ${log_sample_id}"
                    Nextflow.exit(1)
                }

                if (!v && info_field_enum !== SampleMeta.InfoField.LONGITUDINAL_SAMPLE) {
                    log.error "got empty value for info field(${info_field_enum}) for ${log_sample_id}"
                    Nextflow.exit(1)
                }

                info_data[info_field_enum] = v
            }

        return info_data
    }

    private static getFileObject(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    private static void setSampleIds(meta_sample, entry, info_data, log_group_id, log) {

        if (info_data.containsKey(SampleMeta.InfoField.LONGITUDINAL_SAMPLE)) {

            if (meta_sample.containsKey('longitudinal_sample_id') && meta_sample.longitudinal_sample_id != entry.sample_id) {
                log.error "got multiple longitudinal samples for ${log_group_id}: ${entry.sample_id}"
                Nextflow.exit(1)
            }

            meta_sample.longitudinal_sample_id = entry.sample_id

        } else if (meta_sample.containsKey('sample_id') && meta_sample.sample_id != entry.sample_id) {

            log.error "got unexpected sample name for ${log_group_id}: ${entry.sample_id}"
            Nextflow.exit(1)

        } else {

            meta_sample.sample_id = entry.sample_id

        }
    }

    private static void setFastqPaths(meta_sample, entry, info_data, log_sample_id, log) {

        if (!info_data.containsKey(SampleMeta.InfoField.LIBRARY_ID)) {
            log.error "missing 'library_id' info field for ${log_sample_id}"
            Nextflow.exit(1)
        }

        if (!info_data.containsKey(SampleMeta.InfoField.LANE)) {
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
        def fastq_key = [info_data[SampleMeta.InfoField.LIBRARY_ID], info_data[SampleMeta.InfoField.LANE]]

        if (!meta_sample.containsKey(SampleMeta.FileType.FASTQ)) {
            meta_sample[SampleMeta.FileType.FASTQ] = [:]
        }

        if (meta_sample[SampleMeta.FileType.FASTQ].containsKey(fastq_key)) {
            log.error "got duplicate lane + library_id data for ${log_sample_id}: ${fastq_key}"
            Nextflow.exit(1)
        }

        meta_sample[SampleMeta.FileType.FASTQ][fastq_key] = ['fwd': getFileObject(fwd), 'rev': getFileObject(rev)]

    }

    private static void setCramPaths(meta_sample) {

        // CRAMs are passed to hmftools as if they were BAMs, e.g. `-bam_file /path/to/tumor.cram`
        // We therefore set the BAM/BAI path to be the CRAM/CRAI path

        if (meta_sample.containsKey(SampleMeta.FileType.CRAM_REDUX)) {
            meta_sample[SampleMeta.FileType.BAM_REDUX] = meta_sample.remove(SampleMeta.FileType.CRAM_REDUX)
        }

        if (meta_sample.containsKey(SampleMeta.FileType.CRAM)) {
            meta_sample[SampleMeta.FileType.BAM] = meta_sample.remove(SampleMeta.FileType.CRAM)
        }

        // The BAI key is used to store the index for both regular/REDUX CRAMs/BAMs
        if (meta_sample.containsKey(SampleMeta.FileType.CRAI)) {
            meta_sample[SampleMeta.FileType.BAI] = meta_sample.remove(SampleMeta.FileType.CRAI)
        }

    }

    private static void checkAndSetFileIndexes(meta_sample, log) {

        // NOTE(LN): Cast keys to list to avoid ConcurrentModificationException
        meta_sample.keySet().toList().each { key ->

            // NOTE(SW): I was going to use two maps but was unable to get an enum map to compile

            def index_enum
            def index_extension

            if (key === SampleMeta.FileType.BAM || key === SampleMeta.FileType.BAM_REDUX) {
                index_enum = SampleMeta.FileType.BAI
                index_extension = 'bai'
            } else if (key === SampleMeta.FileType.CRAM || key === SampleMeta.FileType.CRAM_REDUX) {
                index_enum = SampleMeta.FileType.CRAI
                index_extension = 'crai'
            } else if (key === SampleMeta.FileType.SAGE_VCF) {
                index_enum = SampleMeta.FileType.SAGE_VCF_TBI
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
                log.error "Could not find index(${index_path}) for file(${file_path})"
                Nextflow.exit(1)
            }

            meta_sample[index_enum] = index_path
        }
    }

    private static void checkRawReadDataExists(meta_sample, log_group_id, log) {

        def missing_any_raw_read_data =
            !meta_sample.containsKey(SampleMeta.FileType.BAM) &&
            !meta_sample.containsKey(SampleMeta.FileType.BAM_REDUX) &&
            !meta_sample.containsKey(SampleMeta.FileType.CRAM) &&
            !meta_sample.containsKey(SampleMeta.FileType.CRAM_REDUX) &&
            !meta_sample.containsKey(SampleMeta.FileType.FASTQ)

        if (missing_any_raw_read_data) {

            log.error "no BAM/CRAM nor BAM_REDUX/CRAM_REDUX nor FASTQ files provided for ${log_group_id}\n\n" +
                    "NB: At least one of these files is required as they are the basis to determine input sample type."
            Nextflow.exit(1)
        }
    }

    private static void setReduxTsvDirIfUnset(meta_sample) {

        if (!meta_sample.containsKey(SampleMeta.FileType.BAM_REDUX)) {
            return
        }

        if(meta_sample.containsKey(SampleMeta.FileType.REDUX_TSV_DIR)) {
            return
        }

        def bam_path = meta_sample[SampleMeta.FileType.BAM_REDUX]
        def bam_dir = bam_path.getParent().toUriString()

        meta_sample[SampleMeta.FileType.REDUX_TSV_DIR] = bam_dir
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
        if (Inputs.hasNormalDna(meta) && !Inputs.hasTumorDna(meta)) {
            log.error "found only normal DNA input for ${meta.group_id} but germline only analysis is not supported"
            Nextflow.exit(1)
        }

        // Do not allow donor sample without normal sample
        if (Inputs.hasDonorDna(meta) && !Inputs.hasNormalDna(meta)) {
            log.error "a donor sample but not normal sample was found for ${meta.group_id}\n\n" +
                "Analysis with a donor sample requires a normal sample."
            Nextflow.exit(1)
        }

        // Apply some required restrictions to targeted mode
        if (run_mode === RunModes.Main.TARGETED) {

            // Do not allow donor DNA
            if (Inputs.hasDonorDna(meta)) {
                log.error "targeted mode is not compatible with the donor DNA BAM/CRAM provided for ${meta.group_id}\n\n" +
                    "The targeted workflow supports only tumor and normal DNA BAM/CRAMs (and tumor RNA BAM/CRAMs for TSO500)"
                Nextflow.exit(1)
            }

            // Do not allow only tumor RNA
            if (Inputs.hasTumorRna(meta) && !Inputs.hasTumorDna(meta)) {
                log.error "targeted mode is not compatible with only tumor RNA provided for ${meta.group_id}\n\n" +
                    "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                    "the configured panel."
                Nextflow.exit(1)
            }
        }
    }

    private static void checkLongitudinalSampleInputs(meta, log) {

        // For purity estimation with WISP, require primary normal DNA BAM when an AMBER directory is provided
        def meta_tumor_dna = meta.getOrDefault([SampleMeta.SampleType.TUMOR, SampleMeta.SequenceType.DNA], [:])
        def longitudinal = meta_tumor_dna.containsKey('longitudinal_sample_id')
        def has_amber_dir = meta_tumor_dna.containsKey(SampleMeta.FileType.AMBER_DIR)
        def has_normal_dna_bam = Inputs.hasNormalDnaBam(meta) || Inputs.hasNormalDnaReduxBam(meta)

        if (longitudinal && has_amber_dir && !has_normal_dna_bam) {
            log.error "AMBER input was provided without the required primary normal DNA BAM for ${meta.group_id}"
            Nextflow.exit(1)
        }

    }
}
