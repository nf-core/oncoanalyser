//
// Input parsing helpers for the nf-core/oncoanalyser pipeline
//

include { getReduxDirAlignment     } from './accessors'
include { Case                     } from './records'
include { FastqFile                } from './records'
include { Sample                   } from './records'
include { FileType                 } from './types'
include { InfoField                } from './types'
include { SampleType               } from './types'
include { SequenceType             } from './types'
include { getEnumFromStringOrFail  } from './utils'
include { getFileObject            } from './utils'

def parseInput(input_fp_str, stub_run, log) {

    if (! input_fp_str) {
        log.error "Missing required --input argument"
        exit 1
    }

    // NOTE(SW): using NF .splitCsv channel operator, hence should be easily interchangable with NF syntax

    def input_fp = getFileObject(input_fp_str)
    def input_entries = nextflow.splitter.SplitterEx.splitCsv(input_fp, [header: true])
    def grouped_entries = input_entries.groupBy { it['case_id'] }
    return grouped_entries.collect { case_id, entries ->
        parseCaseEntry(case_id, entries, stub_run, log)
    }
}

def parseCaseEntry(case_id, entries, stub_run, log) {

    def patient_ids = entries.collect { it.patient_id }.unique()
    if (patient_ids.size() > 1) {
        log.error "got multiple subject names for ${case_id}: ${patient_ids}"
        exit 1
    }
    def patient_id = patient_ids[0]

    def ctx = [cancer_type: null]
    def sample_builders = [:]

    entries.each { entry ->
        parseSampleEntry(case_id, patient_id, ctx, entry, sample_builders, log)
    }

    def samples = sample_builders.values().collect { b ->
        promoteAlignmentFiles(b)
        resolveReduxInputs(b, case_id, stub_run, log)
        checkAlignmentIndexes(b, case_id, stub_run, log)
        record(sample_id: b.sample_id, case_id: case_id, patient_id: patient_id, sample_type: b.sample_type, sequence_type: b.sequence_type, files: b.files, generate_redux_tsvs_only: b.generate_redux_tsvs_only)
    }

    // NOTE(SW): dna_rna samples are folded into the DNA lists until the dna_rna overhaul lands
    def normal_dna = samples.findAll { it.sample_type == SampleType.NORMAL }
    def donor_dna = samples.findAll { it.sample_type == SampleType.DONOR }
    def tumor_dna = samples.findAll { it.sample_type == SampleType.TUMOR && it.sequence_type != SequenceType.RNA }
    def tumor_rna = samples.findAll { it.sample_type == SampleType.TUMOR && it.sequence_type == SequenceType.RNA }
    def longitudinal = samples.findAll { it.sample_type == SampleType.LONGITUDINAL }

    return record(case_id: case_id, patient_id: patient_id, cancer_type: ctx.cancer_type, normal_dna_samples: normal_dna, donor_dna_samples: donor_dna, tumor_dna_samples: tumor_dna, tumor_rna_samples: tumor_rna, longitudinal_samples: longitudinal)
}

def parseSampleEntry(case_id, patient_id, ctx, entry, sample_builders, log) {

    // Sample type
    def sample_type_enum = getEnumFromStringOrFail(entry.sample_type, SampleType, 'sample type', log)

    // Sequence type
    def sequence_type_enum = getEnumFromStringOrFail(entry.sequence_type, SequenceType, 'sequence type', log)

    // Filetype
    def filetype_enum = getEnumFromStringOrFail(entry.filetype, FileType, 'file type', log)

    def sample_key = "${sample_type_enum}_${sequence_type_enum}_${entry.sample_id}"
    def b = sample_builders[sample_key]
    if (! b) {
        b = [sample_id: entry.sample_id, sample_type: sample_type_enum, sequence_type: sequence_type_enum, files: [:], generate_redux_tsvs_only: false]
        sample_builders[sample_key] = b
    }

    // Info data
    def info_data = [:]
    if (entry.containsKey('info') && entry.info) {
        info_data = parseInfoFields(entry.info, case_id, sample_type_enum, sequence_type_enum, log)

        if (info_data.containsKey(InfoField.CANCER_TYPE)) {
            ctx.cancer_type = info_data[InfoField.CANCER_TYPE]
        }

        if (info_data.containsKey(InfoField.GENERATE_REDUX_TSVS_ONLY)) {
            b.generate_redux_tsvs_only = true
        }

        // Only allow READ_GROUP_OVERRIDES for FASTQ
        if (filetype_enum != FileType.FASTQ && info_data.containsKey(InfoField.READ_GROUP_OVERRIDES)) {
            log.error "The read_group info field is only applicable to FASTQ input but got '${entry.filetype}' for ${case_id} ${sample_type_enum}/${sequence_type_enum}"
            exit 1
        }
    }

    // Handle inputs appropriately
    if (filetype_enum == FileType.FASTQ) {
        parseFastqFile(b, entry, info_data, case_id, log)
    } else {
        if (b.files.containsKey(filetype_enum)) {
            log.error "got duplicate file for ${case_id} ${sample_type_enum}/${sequence_type_enum} ${b.sample_id}: ${filetype_enum}"
            exit 1
        }
        b.files[filetype_enum] = getFileObject(entry.filepath)
    }
}

def parseFastqFile(b, entry, info_data, case_id, log) {

    if (! info_data.containsKey(InfoField.LIBRARY_ID)) {
        log.error "missing 'library_id' info field for ${case_id} ${b.sample_type}/${b.sequence_type} ${b.sample_id}"
        exit 1
    }

    if (! info_data.containsKey(InfoField.LANE)) {
        log.error "missing 'lane' info field for ${case_id} ${b.sample_type}/${b.sequence_type} ${b.sample_id}"
        exit 1
    }

    def fastq_entries = entry.filepath.tokenize(';')

    if (fastq_entries.size() != 1 && fastq_entries.size() != 2) {
        log.error "expected 1 (single-end) or 2 (paired-end) FASTQ files delimited by ';' (i.e. '<fwd>' or '<fwd>;<rev>') but found ${fastq_entries.size} " +
            " files for ${case_id} ${b.sample_type}/${b.sequence_type} ${b.sample_id}"
        exit 1
    }

    def single_end = fastq_entries.size() == 1
    def fwd = getFileObject(fastq_entries[0])
    def rev = single_end ? null : getFileObject(fastq_entries[1])

    def rg_fields = [:]
    if (info_data.containsKey(InfoField.READ_GROUP_OVERRIDES)) {
        rg_fields = parse_read_group_info(info_data[InfoField.READ_GROUP_OVERRIDES], log)
    }

    if (! b.files.containsKey(FileType.FASTQ)) {
        b.files[FileType.FASTQ] = []
    }

    def fastq = record(
        read_fwd: fwd,
        read_rev: rev,
        single_end: single_end,
        library_id: info_data[InfoField.LIBRARY_ID],
        lane: info_data[InfoField.LANE],
        flowcell: info_data.getOrDefault(InfoField.FLOWCELL, null),
        rg_fields: rg_fields,
    )

    def duplicate = b.files[FileType.FASTQ].find { existing ->
        existing.library_id == fastq.library_id && existing.lane == fastq.lane && existing.flowcell == fastq.flowcell
    }
    if (duplicate) {
        log.error "got duplicate lane + library_id data for ${case_id} ${b.sample_type}/${b.sequence_type} ${b.sample_id}"
        exit 1
    }

    b.files[FileType.FASTQ] << fastq
}

def parseInfoFields(info_str, case_id, sample_type_enum, sequence_type_enum, log) {
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

            def info_field_enum = getEnumFromStringOrFail(k, InfoField, 'info field', log)

            if (info_data.containsKey(info_field_enum)) {
                log.error "got duplicate info field for ${case_id} ${sample_type_enum}/${sequence_type_enum}: ${info_field_enum}"
                exit 1
            }

            if (! v && info_field_enum != InfoField.LONGITUDINAL_SAMPLE && info_field_enum != InfoField.GENERATE_REDUX_TSVS_ONLY) {
                log.error "got empty value for ${case_id} ${sample_type_enum}/${sequence_type_enum} ${info_field_enum}"
                exit 1
            }

            info_data[info_field_enum] = v
        }

    return info_data
}

def promoteAlignmentFiles(b) {
    def files = b.files

    if (files.containsKey(FileType.BAM)) {
        files[FileType.ALN] = files.remove(FileType.BAM)
    }

    if (files.containsKey(FileType.CRAM)) {
        files[FileType.ALN] = files.remove(FileType.CRAM)
    }

    if (files.containsKey(FileType.BAM_REDUX)) {
        files[FileType.ALN_REDUX] = files.remove(FileType.BAM_REDUX)
    }

    if (files.containsKey(FileType.CRAM_REDUX)) {
        files[FileType.ALN_REDUX] = files.remove(FileType.CRAM_REDUX)
    }

    if (files.containsKey(FileType.BAI)) {
        files[FileType.IDX] = files.remove(FileType.BAI)
    }

    if (files.containsKey(FileType.CRAI)) {
        files[FileType.IDX] = files.remove(FileType.CRAI)
    }
}

def resolveReduxInputs(b, case_id, stub_run, log) {

    if (stub_run) {
        return
    }

    def files = b.files
    def sample_id = b.sample_id

    if (files.containsKey(FileType.ALN_REDUX) && files.containsKey(FileType.REDUX_DIR)) {
        log.error "expected either REDUX directory or REDUX alignment but got both for ${case_id} ${sample_id}"
        exit 1
    }

    def redux_input
    def redux_aln
    def redux_dir
    if (files.containsKey(FileType.ALN_REDUX)) {

        redux_aln = files[FileType.ALN_REDUX]
        redux_input = redux_aln
        redux_dir = redux_aln.parent
        if (! redux_input.isFile()) {
            log.error "didn't receive file as REDUX alignment for ${case_id} ${sample_id}: ${redux_input}"
            exit 1
        }

    } else if (files.containsKey(FileType.REDUX_DIR)) {

        redux_dir = files[FileType.REDUX_DIR]
        redux_input = redux_dir
        redux_aln = getReduxDirAlignment(sample_id, redux_dir)[0]
        if (! redux_input.isDirectory()) {
            log.error "didn't receive directory as REDUX directory for ${case_id} ${sample_id}: ${redux_input}"
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

    def generate_tsvs_only = b.generate_redux_tsvs_only

    if (! has_redux_tsvs && ! generate_tsvs_only) {

        log.error "no REDUX TSVs provided or found for ${case_id} ${sample_id}: ${redux_input}"
        exit 1

    }

    if (! has_colocated_index) {

        if (files.containsKey(FileType.REDUX_DIR)) {
            log.error "required index not located in REDUX directory: ${case_id} ${sample_id}: ${redux_input}"
            exit 1
        }

    }

    if (has_colocated_index && ! files.containsKey(FileType.REDUX_DIR)) {
        files.remove(FileType.ALN_REDUX)
        files[FileType.REDUX_DIR] = redux_dir
    }

    if (! has_redux_tsvs && generate_tsvs_only) {
        files.remove(FileType.REDUX_DIR)
        files[FileType.ALN_REDUX] = redux_aln
    }

    if (files.containsKey(FileType.ALN_REDUX) && files.containsKey(FileType.IDX) && ! has_redux_tsvs && ! generate_tsvs_only) {
        log.error "REDUX alignments without colocated TSVs requires generate_redux_tsvs_only to be set in the samplesheet: ${case_id} ${sample_id}: ${redux_input}"
        exit 1
    }

    if (files.containsKey(FileType.REDUX_DIR) && has_redux_tsvs && generate_tsvs_only) {
        log.warn "REDUX directory already contains TSVs, ignoring generate_redux_tsvs_only flag for: ${case_id} ${sample_id}: ${redux_input}"
    }
}

def checkAlignmentIndexes(b, case_id, stub_run, log) {

    if (stub_run) {
        return
    }

    def files = b.files
    def sample_id = b.sample_id

    files.keySet().toList().each { key ->

        def aln
        if (key == FileType.ALN || key == FileType.ALN_REDUX) {
            aln = files[key]
        } else if (key == FileType.REDUX_DIR) {
            def d = getReduxDirAlignment(sample_id, files[key])
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
            log.error "got bad alignment type for ${case_id} ${b.sample_type}/${b.sequence_type}: ${key}: ${aln}"
            exit 1
        }

        if (files.containsKey(FileType.IDX) && files[FileType.IDX].name.endsWith(index_ext)) {
            return
        }

        def fp = aln.toUriString()
        def index_fp = nextflow.Nextflow.file("${fp}.${index_ext}")

        if (! index_fp.exists() && ! stub_run) {
            log.error "no index provided or found for ${case_id} ${b.sample_type}/${b.sequence_type}: ${key}: ${fp}"
            exit 1
        }

        if (key == FileType.ALN || key == FileType.ALN_REDUX) {
            files[FileType.IDX] = index_fp
        }
    }
}

def parse_read_group_info(rg_info_raw, log) {
    def escape_char = "\u0000"
    def validate_rg_tags = ['BC', 'CN', 'DS', 'DT', 'FO', 'ID', 'KS', 'LB', 'PG', 'PI', 'PL', 'PM', 'PU', 'SM']

    def fields = [:]
    def rg_info_escaped = rg_info_raw.replace('||', escape_char)
    rg_info_escaped.split('\\|').each { field_str_escaped ->
        def field_str = field_str_escaped.replace(escape_char, '|')
        if (! field_str.contains('=')) {
            log.error "Received bad read group field (must be in format `<name>=<value>`): ${field_str}"
            exit 1
        }

        def (name, value) = field_str.split('=', 2)
        if (! validate_rg_tags.contains(name)) {
            log.error "Received bad read group tag '${name}' in: ${rg_info_raw}"
            exit 1
        }

        if (! value) {
            log.error "Received empty read group value for '${name}' in: ${rg_info_raw}"
            exit 1
        }

        fields[name] = value
    }

    return fields
}
