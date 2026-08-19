//
// Input parsing helpers for the nf-core/oncoanalyser pipeline
//

include { getEnumFromStringOrFail; getFileObject; getReduxDirAlignment; hasNormalDnaBam; hasNormalDnaReduxInput; parse_read_group_info } from './utils'

def parseInput(input_fp_str, stub_run, log) {

    if (! input_fp_str) {
        log.error "Missing required --input argument"
        exit 1
    }

    // NOTE(SW): using NF .splitCsv channel operator, hence should be easily interchangable with NF syntax

    def input_fp = getFileObject(input_fp_str)
    def input_entries = nextflow.splitter.SplitterEx.splitCsv(input_fp, [header: true])
    def grouped_entries = input_entries.groupBy { it['group_id'] }
    def inputs = grouped_entries.collect { group_id, entries ->
        parseCaseEntry(group_id, entries, stub_run, log)
    }

    return inputs
}

def parseCaseEntry(group_id, entries, stub_run, log) {

    def meta = [group_id: group_id]
    def sample_keys = [] as Set

    // Process each entry
    entries.each {
        parseSampleEntry(meta, sample_keys, it, log)
    }

    // Set BAM / CRAM as generic alignments, along with their index
    promoteAlignmentFiles(meta, sample_keys)

    // Handle REDUX alignments; require TSVs to be present, and promote to directory when files for this sample are found.
    // Allow multiple samples to be in the same dir (e.g. both tumor and normal sample).
    // Existing helpers (getReduxDirAlignment, getReduxTsvs) already resolve files by exact {sample_id}.redux.* filename.
    resolveReduxInputs(meta, sample_keys, stub_run, log)

    // Check that required indexes are provided or are accessible
    checkAlignmentIndexes(meta, sample_keys, stub_run, log)

    // For purity estimation with WISP, require primary normal DNA BAM when an AMBER directory is provided
    def meta_tumor_dna = meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
    def longitudinal = meta_tumor_dna.containsKey('longitudinal_sample_id')
    def has_amber_dir = meta_tumor_dna.containsKey(Constants.FileType.AMBER_DIR)
    def has_normal_dna_bam = hasNormalDnaBam(meta) || hasNormalDnaReduxInput(meta)

    if (longitudinal && has_amber_dir && ! has_normal_dna_bam) {
        log.error "AMBER input was provided without the required primary normal DNA BAM for ${meta.group_id}"
        exit 1
    }

    return meta
}

def parseSampleEntry(meta, sample_keys, entry, log) {

    // Add subject id if absent or check if current matches existing
    if (meta.containsKey('subject_id') && meta.subject_id != entry.subject_id) {
        log.error "got unexpected subject name for ${meta.group_id} ${meta.subject_id}: ${entry.subject_id}"
        exit 1
    } else {
        meta.subject_id = entry.subject_id
    }

    // Sample type
    def sample_type_enum = getEnumFromStringOrFail(entry.sample_type, Constants.SampleType, 'sample type', log)

    // Sequence type
    def sequence_type_enum = getEnumFromStringOrFail(entry.sequence_type, Constants.SequenceType, 'sequence type', log)

    // Filetype
    def filetype_enum = getEnumFromStringOrFail(entry.filetype, Constants.FileType, 'file type', log)

    def sample_key = [sample_type_enum, sequence_type_enum]
    def meta_sample = meta.get(sample_key, [:])

    // Info data
    def info_data = [:]
    if (entry.containsKey('info')) {
        // Parse
        info_data = parseInfoFields(entry.info, meta.group_id, sample_type_enum, sequence_type_enum, log)

        // Process
        if (info_data.containsKey(Constants.InfoField.CANCER_TYPE)) {
            meta[Constants.InfoField.CANCER_TYPE] = info_data[Constants.InfoField.CANCER_TYPE]
        }

        if (info_data.containsKey(Constants.InfoField.GENERATE_REDUX_TSVS_ONLY)) {
            meta_sample[Constants.InfoField.GENERATE_REDUX_TSVS_ONLY] = true
        }

        // Only allow READ_GROUP_OVERRIDES for FASTQ
        if (filetype_enum != Constants.FileType.FASTQ && info_data.containsKey(Constants.InfoField.READ_GROUP_OVERRIDES)) {
            log.error "The read_group info field is only applicable to FASTQ input but got '${entry.filetype}' for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}"
            exit 1
        }

    }

    if (info_data.containsKey(Constants.InfoField.LONGITUDINAL_SAMPLE)) {

        if (meta_sample.containsKey('longitudinal_sample_id') && meta_sample.longitudinal_sample_id != entry.sample_id) {
            log.error "got multiple longitudinal samples for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}: ${entry.sample_id}"
            exit 1
        }

        meta_sample.longitudinal_sample_id = entry.sample_id

    } else if (meta_sample.containsKey('sample_id') && meta_sample.sample_id != entry.sample_id) {

        log.error "got unexpected sample name for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}: ${entry.sample_id}"
        exit 1

    } else {

        meta_sample.sample_id = entry.sample_id

    }

    // Disallow AMBER, COBALT, SAGE_APPEND inputs for longitudinal samples; these would clash with primary inputs
    if (info_data.containsKey(Constants.InfoField.LONGITUDINAL_SAMPLE)) {
        def longitudinal_disallowed_input_list = [Constants.FileType.AMBER_DIR, Constants.FileType.COBALT_DIR, Constants.FileType.SAGE_APPEND_DIR]
        def disallowed_inputs = longitudinal_disallowed_input_list.findAll { e -> e == filetype_enum }
        if (disallowed_inputs) {
            log.error "got disallowed ${filetype_enum} input for longitudinal sample ${meta.group_id} ${meta_sample.sample_id} ${sample_type_enum}/${sequence_type_enum}"
            exit 1
        }
    }

    // Filetype uniqueness
    if (meta_sample.containsKey(filetype_enum) && filetype_enum != Constants.FileType.FASTQ) {
        log.error "got duplicate file for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}: ${filetype_enum}"
        exit 1
    }

    // Handle inputs appropriately
    if (filetype_enum == Constants.FileType.FASTQ) {

        if (! info_data.containsKey(Constants.InfoField.LIBRARY_ID)) {
            log.error "missing 'library_id' info field for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}"
            exit 1
        }

        if (! info_data.containsKey(Constants.InfoField.LANE)) {
            log.error "missing 'lane' info field for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}"
            exit 1
        }

        def fastq_entries = entry.filepath.tokenize(';')

        if (fastq_entries.size() != 2) {
            log.error "expected exactly 2 FASTQ files delimited by ';' (i.e. '<fwd>;<rev>') but found ${fastq_entries.size} " +
                " files for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}"
            exit 1
        }

        def (fwd, rev) = fastq_entries
        def fastq_key = [info_data[Constants.InfoField.LIBRARY_ID], info_data[Constants.InfoField.LANE], info_data.getOrDefault(Constants.InfoField.FLOWCELL, null)]

        if (! meta_sample.containsKey(filetype_enum)) {
            meta_sample[filetype_enum] = [:]
        }

        if (meta_sample[filetype_enum].containsKey(fastq_key)) {
            log.error "got duplicate lane + library_id data for ${meta.group_id} ${sample_type_enum}/${sequence_type_enum}: ${fastq_key}"
            exit 1
        }

        def rg_fields = [:]
        if (info_data.containsKey(Constants.InfoField.READ_GROUP_OVERRIDES)) {
            rg_fields = parse_read_group_info(info_data[Constants.InfoField.READ_GROUP_OVERRIDES], log)
        }

        meta_sample[filetype_enum][fastq_key] = ['fwd': getFileObject(fwd), 'rev': getFileObject(rev), 'rg_fields': rg_fields]

    } else {

        meta_sample[filetype_enum] = getFileObject(entry.filepath)

    }

    // Record sample key to simplify iteration later on
    sample_keys << sample_key
}

def parseInfoFields(info_str, group_id, sample_type_enum, sequence_type_enum, log) {
    def info_data = [:]

    info_str
        .tokenize(';')
        .each { e ->

            def k
            def v
            if (e.contains(':')) {
               (k, v) = e.split(':', 2)
            } else {
               k = e
            }

            def info_field_enum = getEnumFromStringOrFail(k, Constants.InfoField, 'info field', log)

            if (info_data.containsKey(info_field_enum)) {
                log.error "got duplicate info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${info_field_enum}"
                exit 1
            }

            if (! v && info_field_enum != Constants.InfoField.LONGITUDINAL_SAMPLE && info_field_enum != Constants.InfoField.GENERATE_REDUX_TSVS_ONLY) {
                log.error "got empty value for ${group_id} ${sample_type_enum}/${sequence_type_enum} ${info_field_enum}"
                exit 1
            }

            info_data[info_field_enum] = v
        }

    return info_data
}

def promoteAlignmentFiles(meta, sample_keys) {
    sample_keys.each { sample_key ->

        def meta_sample = meta[sample_key]

        if (meta_sample.containsKey(Constants.FileType.BAM)) {
            meta_sample[Constants.FileType.ALN] = meta_sample.remove(Constants.FileType.BAM)
        }

        if (meta_sample.containsKey(Constants.FileType.CRAM)) {
            meta_sample[Constants.FileType.ALN] = meta_sample.remove(Constants.FileType.CRAM)
        }

        if (meta_sample.containsKey(Constants.FileType.BAM_REDUX)) {
            meta_sample[Constants.FileType.ALN_REDUX] = meta_sample.remove(Constants.FileType.BAM_REDUX)
        }

        if (meta_sample.containsKey(Constants.FileType.CRAM_REDUX)) {
            meta_sample[Constants.FileType.ALN_REDUX] = meta_sample.remove(Constants.FileType.CRAM_REDUX)
        }

        if (meta_sample.containsKey(Constants.FileType.BAI)) {
            meta_sample[Constants.FileType.IDX] = meta_sample.remove(Constants.FileType.BAI)
        }

        if (meta_sample.containsKey(Constants.FileType.CRAI)) {
            meta_sample[Constants.FileType.IDX] = meta_sample.remove(Constants.FileType.CRAI)
        }

    }
}

def resolveReduxInputs(meta, sample_keys, stub_run, log) {
    sample_keys.each { sample_key ->

        def meta_sample = meta[sample_key]
        def is_primary = ! meta_sample.containsKey('longitudinal_sample_id')
        def sample_id = is_primary ? meta_sample.sample_id : meta_sample.longitudinal_sample_id

        if (stub_run) {
            return
        }

        if (meta_sample.containsKey(Constants.FileType.ALN_REDUX) && meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {
            log.error "expected either REDUX directory or REDUX alignment but got both for ${meta.group_id} ${sample_id}"
            exit 1
        }

        def redux_input
        def redux_aln
        def redux_dir
        if (meta_sample.containsKey(Constants.FileType.ALN_REDUX)) {

            redux_aln = meta_sample[Constants.FileType.ALN_REDUX]
            redux_input = redux_aln
            redux_dir = redux_aln.parent
            if (! redux_input.isFile()) {
                log.error "didn't receive file as REDUX alignment for ${meta.group_id} ${sample_id}: ${redux_input}"
                exit 1
            }

        } else if (meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {

            redux_dir = meta_sample[Constants.FileType.REDUX_DIR]
            redux_input = redux_dir
            redux_aln = getReduxDirAlignment(sample_id, redux_dir)[0]
            if (! redux_input.isDirectory()) {
                log.error "didn't receive directory as REDUX directory for ${meta.group_id} ${sample_id}: ${redux_input}"
                exit 1
            }

        } else {

            return

        }

        def has_bqr_tsv = redux_dir.resolve("${sample_id}.redux.bqr.tsv").exists()
        def has_jitter_tsv = redux_dir.resolve("${sample_id}.redux.jitter_params.tsv").exists()
        def has_ms_tsv = redux_dir.resolve("${sample_id}.redux.ms_table.tsv.gz").exists()
        def has_redux_tsvs = has_bqr_tsv && has_jitter_tsv && has_ms_tsv

        def has_colocated_index = nextflow.Nextflow.file("${redux_aln.toUriString()}.bai").exists() || nextflow.Nextflow.file("${redux_aln.toUriString()}.crai").exists()

        def generate_tsvs_only = meta_sample.getOrDefault(Constants.InfoField.GENERATE_REDUX_TSVS_ONLY, false)

        if (! has_redux_tsvs && ! generate_tsvs_only) {

            log.error "no REDUX TSVs provided or found for ${meta.group_id} ${sample_id}: ${redux_input}"
            exit 1

        }

        if (! has_colocated_index) {

            if (meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {
                log.error "required index not located in REDUX directory: ${meta.group_id} ${sample_id}: ${redux_input}"
                exit 1
            }

        }

        if (has_colocated_index && ! meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {
            meta_sample.remove(Constants.FileType.ALN_REDUX)
            meta_sample[Constants.FileType.REDUX_DIR] = redux_dir
        }

        if (! has_redux_tsvs && generate_tsvs_only) {
            meta_sample.remove(Constants.FileType.REDUX_DIR)
            meta_sample[Constants.FileType.ALN_REDUX] = redux_aln
        }

        if (meta_sample.containsKey(Constants.FileType.ALN_REDUX) && meta_sample.containsKey(Constants.FileType.IDX) && ! has_redux_tsvs && ! generate_tsvs_only) {
            log.error "REDUX alignments without colocated TSVs requires generate_redux_tsvs_only to be set in the samplesheet: ${meta.group_id} ${sample_id}: ${redux_input}"
            exit 1
        }

        if (meta_sample.containsKey(Constants.FileType.REDUX_DIR) && has_redux_tsvs && generate_tsvs_only) {
            log.warn "REDUX directory already contains TSVs, ignoring generate_redux_tsvs_only flag for: ${meta.group_id} ${sample_id}: ${redux_input}"
        }

    }
}

def checkAlignmentIndexes(meta, sample_keys, stub_run, log) {
    sample_keys.each { sample_key ->

        meta[sample_key]*.key.each { key ->

            def meta_sample = meta[sample_key]
            def is_primary = ! meta_sample.containsKey('longitudinal_sample_id')
            def sample_id = is_primary ? meta_sample.sample_id : meta_sample.longitudinal_sample_id

            if (stub_run) {
                return
            }

            def aln
            if (key == Constants.FileType.ALN || key == Constants.FileType.ALN_REDUX) {
                aln = meta_sample[key]
            } else if (key == Constants.FileType.REDUX_DIR) {
                def d = getReduxDirAlignment(sample_id, meta_sample[key])
                aln = d[0]
            } else {
                return
            }

            def index_ext
            if (aln.name.endsWith('.bam')) {
                index_ext = 'bai'
            } else if (aln.name.endsWith('.cram')) {
                index_ext = 'crai'
            } else {
                def (sample_type, sequence_type) = sample_key
                log.error "got bad alignment type for ${meta.group_id} ${sample_type}/${sequence_type}: ${key}: ${aln}"
                exit 1
            }

            if (meta_sample.containsKey(Constants.FileType.IDX) && meta_sample[Constants.FileType.IDX].name.endsWith(index_ext)) {
                return
            }

            def fp = aln.toUriString()
            def index_fp = nextflow.Nextflow.file("${fp}.${index_ext}")

            if (! index_fp.exists() && ! stub_run) {
                def (sample_type, sequence_type) = sample_key
                log.error "no index provided or found for ${meta.group_id} ${sample_type}/${sequence_type}: ${key}: ${fp}"
                exit 1
            }

            if (key == Constants.FileType.ALN || key == Constants.FileType.ALN_REDUX) {
                meta_sample[Constants.FileType.IDX] = index_fp
            }
        }
    }
}
