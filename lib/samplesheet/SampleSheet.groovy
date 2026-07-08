package samplesheet

import pipeline.PipelineMode
import sample.Inputs
import util.Messages

class SampleSheet {

    public static List<Map> parse(String sample_sheet_path, PipelineMode pipeline_mode, boolean stub_run) {

        if (!sample_sheet_path) {
            throw new IllegalStateException("Missing required --input argument")
        }

        // NOTE(SW): using NF .splitCsv channel operator, hence should be easily interchangable with NF syntax

        def sample_sheet_file = getFileObject(sample_sheet_path)
        def inputs = nextflow.splitter.SplitterEx.splitCsv(sample_sheet_file, [header: true])
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

                entries.each { entry -> createOrUpdateSampleMeta(entry, meta, sample_keys) }

                // Per sample checks once meta_sample objects are fully constructed
                sample_keys.each { sample_key ->

                    def meta_sample = meta[sample_key]

                    // NOTE(LN): The order in which these methods are executed is important
                    checkAndSetFileIndexes(meta_sample, stub_run)
                    setCramPaths(meta_sample)
                    checkRawReadDataExists(meta_sample, group_id)
                    setReduxTsvDirIfUnset(meta_sample)

                }

                // Per group checks
                disallowDuplicateSampleIdsWithinSampleGroup(meta)
                disallowInvalidSampleCombinations(meta, pipeline_mode)
                checkLongitudinalSampleInputs(meta)

                return meta
            }

        return inputs
    }

    private static void createOrUpdateSampleMeta(Map<String, String> entry, Map meta, Set<List> sample_keys) {

        def sample_type = SampleType.fromString(entry.sample_type)
        def sequence_type = SequenceType.fromString(entry.sequence_type)
        def file_type = FileType.fromString(entry.filetype)

        setAndCheckSubjectId(meta, entry)

        def sample_key = [sample_type, sequence_type]
        sample_keys.add(sample_key) // Store to iterate over later

        def meta_sample = meta.get(sample_key, [:])

        // Handle info field
        def info_data = parseInfoField(entry)

        if (info_data.containsKey(InfoField.CANCER_TYPE)) {
            meta[InfoField.CANCER_TYPE] = info_data[InfoField.CANCER_TYPE]
        }

        setSampleIds(meta_sample, entry, info_data)

        // Set file paths
        if (meta_sample.containsKey(file_type) & file_type != FileType.FASTQ) {
            throw new IllegalStateException("Got duplicate filetype(${file_type}) for group_id(${entry.group_id}) sample_id(${entry.sample_id})")
        }

        if (file_type === FileType.FASTQ) {
            setFastqPaths(meta_sample, entry, info_data)
        } else {
            meta_sample[file_type] = getFileObject(entry.filepath)
        }
    }

    private static void setAndCheckSubjectId(meta, entry) {

        if(!meta.containsKey('subject_id')) {
            meta.subject_id = entry.subject_id
        }

        if (meta.subject_id != entry.subject_id) {
            throw new IllegalStateException(
                "Got unexpected subject_id(${entry.subject_id}) for group_id(${entry.group_id}). " +
                "All samples within a group must have the same subject_id"
            )
        }
    }

    private static Map parseInfoField(Map<String, String> entry) {
        def info_data = [:]

        if (!entry.containsKey('info')) {
            return info_data
        }

        entry.info
            .tokenize(';')
            .each { info_item ->
                def (info_key, info_value) = info_item.tokenize(':')
                def info_field = InfoField.fromString(info_key)

                if (info_data.containsKey(info_field)) {
                    throw new IllegalStateException("Got duplicate info field(${info_field}) for group_id(${entry.group_id}) sample_id(${entry.sample_id})")
                }

                if (!info_value && info_field !== InfoField.LONGITUDINAL_SAMPLE) {
                    throw new IllegalStateException("Got empty value for info field(${info_field}) for group_id(${entry.group_id}) sample_id(${entry.sample_id})")
                }

                info_data[info_field] = info_value
            }

        return info_data
    }

    private static getFileObject(String path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    private static void setSampleIds(Map meta_sample, Map<String, String> entry, Map<InfoField, String> info_data) {

        if (info_data.containsKey(InfoField.LONGITUDINAL_SAMPLE)) {

            if (meta_sample.containsKey('longitudinal_sample_id') && meta_sample.longitudinal_sample_id != entry.sample_id) {
                throw new IllegalStateException("Got multiple longitudinal samples for group_id(${entry.group_id}) - this is currently unsupported")
            }

            meta_sample.longitudinal_sample_id = entry.sample_id

        } else if (meta_sample.containsKey('sample_id') && meta_sample.sample_id != entry.sample_id) {

            throw new IllegalStateException("Got unexpected sample_id(${entry.sample_id}) for group_id(${entry.group_id})")

        } else {

            meta_sample.sample_id = entry.sample_id

        }
    }

    private static void setFastqPaths(Map meta_sample, Map<String, String> entry, Map<InfoField, String> info_data) {

        if (!info_data.containsKey(InfoField.LIBRARY_ID)) {
            throw new IllegalStateException("Missing info field(library_id) for group_id(${entry.group_id}) sample_id(${entry.sample_id})")
        }

        if (!info_data.containsKey(InfoField.LANE)) {
            throw new IllegalStateException("Missing 'lane' info field for group_id(${entry.group_id}) sample_id(${entry.sample_id})")
        }

        def fastq_entries = entry.filepath.tokenize(';')

        if (fastq_entries.size() != 2) {
            throw new IllegalStateException(
                "Expected exactly 2 FASTQ files delimited by ';' (i.e. '<fwd>;<rev>') but found " +
                "${fastq_entries.size} files for group_id(${entry.group_id}) sample_id(${entry.sample_id})"
            )
        }

        def (fwd, rev) = fastq_entries
        def fastq_key = [info_data[InfoField.LIBRARY_ID], info_data[InfoField.LANE]]

        if (!meta_sample.containsKey(FileType.FASTQ)) {
            meta_sample[FileType.FASTQ] = [:]
        }

        if (meta_sample[FileType.FASTQ].containsKey(fastq_key)) {
            throw new IllegalStateException("Got duplicate lane + library_id data for group_id(${entry.group_id}) sample_id(${entry.sample_id}): ${fastq_key}")
        }

        meta_sample[FileType.FASTQ][fastq_key] = ['fwd': getFileObject(fwd), 'rev': getFileObject(rev)]

    }

    private static void setCramPaths(Map meta_sample) {

        // CRAMs are passed to hmftools as if they were BAMs, e.g. `-bam_file /path/to/tumor.cram`
        // We therefore set the BAM/BAI path to be the CRAM/CRAI path

        if (meta_sample.containsKey(FileType.CRAM_REDUX)) {
            meta_sample[FileType.BAM_REDUX] = meta_sample.remove(FileType.CRAM_REDUX)
        }

        if (meta_sample.containsKey(FileType.CRAM)) {
            meta_sample[FileType.BAM] = meta_sample.remove(FileType.CRAM)
        }

        // The BAI key is used to store the index for both regular/REDUX CRAMs/BAMs
        if (meta_sample.containsKey(FileType.CRAI)) {
            meta_sample[FileType.BAI] = meta_sample.remove(FileType.CRAI)
        }

    }

    private static void checkAndSetFileIndexes(Map meta_sample, boolean stub_run) {

        // NOTE(LN): Cast keys to list to avoid ConcurrentModificationException
        meta_sample.keySet().toList().each { key ->

            // NOTE(SW): I was going to use two maps but was unable to get an enum map to compile

            def index_enum
            def index_extension

            if (key === FileType.BAM || key === FileType.BAM_REDUX) {
                index_enum = FileType.BAI
                index_extension = 'bai'
            } else if (key === FileType.CRAM || key === FileType.CRAM_REDUX) {
                index_enum = FileType.CRAI
                index_extension = 'crai'
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

            if (!index_path.exists() && !stub_run) {
                throw new IllegalStateException("Could not find index(${index_path}) for file(${file_path})")
            }

            meta_sample[index_enum] = index_path
        }
    }

    private static void checkRawReadDataExists(Map meta_sample, String group_id) {

        def missing_any_raw_read_data =
            !meta_sample.containsKey(FileType.BAM) &&
            !meta_sample.containsKey(FileType.BAM_REDUX) &&
            !meta_sample.containsKey(FileType.CRAM) &&
            !meta_sample.containsKey(FileType.CRAM_REDUX) &&
            !meta_sample.containsKey(FileType.FASTQ)

        if (missing_any_raw_read_data) {
            Messages.error(
                "No BAM/CRAM nor BAM_REDUX/CRAM_REDUX nor FASTQ files provided for group_id(${group_id})",
                "",
                "NB: At least one of these files is required as they are the basis to determine input sample type."
            )
        }
    }

    private static void setReduxTsvDirIfUnset(Map meta_sample) {

        if (!meta_sample.containsKey(FileType.BAM_REDUX)) {
            return
        }

        if(meta_sample.containsKey(FileType.REDUX_TSV_DIR)) {
            return
        }

        def bam_path = meta_sample[FileType.BAM_REDUX]
        def bam_dir = bam_path.getParent().toUriString()

        meta_sample[FileType.REDUX_TSV_DIR] = bam_dir
    }

    private static void disallowDuplicateSampleIdsWithinSampleGroup(Map meta) {

        def sample_ids_seen = [] as Set

        for (maybe_meta_sample in meta.values()) {

            def is_meta_sample = (maybe_meta_sample instanceof Map) && maybe_meta_sample.containsKey('sample_id')
            if(!is_meta_sample)
                continue

            def sample_id = maybe_meta_sample.sample_id

            if(sample_ids_seen.contains(sample_id)) {
                throw new IllegalStateException("Duplicate sample_id(${sample_id}) found for group_id(${meta.group_id})")
            }

            sample_ids_seen.add(sample_id)
        }
    }

    private static void disallowInvalidSampleCombinations(Map meta, PipelineMode pipeline_mode) {

        // Do not allow normal DNA only
        if (Inputs.hasNormalDna(meta) && !Inputs.hasTumorDna(meta)) {
            throw new IllegalStateException("Found only normal DNA input for group_id(${meta.group_id}) but germline only analysis is not supported")
        }

        // Do not allow donor sample without normal sample
        if (Inputs.hasDonorDna(meta) && !Inputs.hasNormalDna(meta)) {
            throw new IllegalStateException("Donor sample provided without normal sample for group_id(${meta.group_id}).")
        }

        // Apply some required restrictions to targeted mode
        if (pipeline_mode === PipelineMode.TARGETED) {

            // Do not allow donor DNA
            if (Inputs.hasDonorDna(meta)) {
                throw new IllegalStateException("Targeted mode is not compatible with the donor DNA BAM/CRAM provided for group_id(${meta.group_id})")
            }

            // Do not allow only tumor RNA
            if (Inputs.hasTumorRna(meta) && !Inputs.hasTumorDna(meta)) {

                Messages.error(
                    "Targeted mode is not compatible with only tumor RNA provided for group_id(${meta.group_id})",
                    "",
                    "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on ",
                    "the configured panel."
                )
            }
        }
    }

    private static void checkLongitudinalSampleInputs(Map meta) {

        // For purity estimation with WISP, require primary normal DNA BAM when an AMBER directory is provided
        def meta_tumor_dna = meta.getOrDefault([SampleType.TUMOR, SequenceType.DNA], [:])
        def longitudinal = meta_tumor_dna.containsKey('longitudinal_sample_id')
        def has_amber_dir = meta_tumor_dna.containsKey(FileType.AMBER_DIR)
        def has_normal_dna_bam = Inputs.hasNormalDnaReduxBam(meta)

        if (longitudinal && has_amber_dir && !has_normal_dna_bam) {
            Messages.error(
                "group_id(${meta.group_id}) - In purity estimate mode, if AMBER dir is provided,",
                "the primary normal DNA REDUX BAM should also be provided"
            )
        }

    }
}
