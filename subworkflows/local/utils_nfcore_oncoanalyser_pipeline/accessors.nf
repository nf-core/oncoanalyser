//
// Sample and file accessors for the nf-core/oncoanalyser pipeline
//

include { FileType     } from './types'
include { SequenceType } from './types'

// All samples in a case, across every sample list
def getSamples(case_record) {
    return case_record.normal_dna_samples + case_record.donor_dna_samples + case_record.tumor_dna_samples + case_record.tumor_rna_samples + case_record.longitudinal_samples
}

// All samples matching the given sample types and sequence types
def getSamples(case_record, sampletypes, sequencetypes) {
    return getSamples(case_record).findAll { s -> sampletypes.contains(s.sample_type) && sequencetypes.contains(s.sequence_type) }
}

def hasAlignmentInput(sample) {
    def files = sample.files
    return files.containsKey(FileType.FASTQ) ||
        files.containsKey(FileType.ALN) ||
        files.containsKey(FileType.ALN_REDUX) ||
        files.containsKey(FileType.REDUX_DIR)
}

// Sample records (singular: first matching sample, or null)
def getTumorDnaSample(case_record) {
    def samples = case_record.tumor_dna_samples.findAll { it.sequence_type == SequenceType.DNA }
    return samples ? samples[0] : null
}

def getTumorRnaSample(case_record) {
    return case_record.tumor_rna_samples ? case_record.tumor_rna_samples[0] : null
}

def getNormalDnaSample(case_record) {
    def samples = case_record.normal_dna_samples.findAll { it.sequence_type == SequenceType.DNA }
    return samples ? samples[0] : null
}

def getDonorDnaSample(case_record) {
    def samples = case_record.donor_dna_samples.findAll { it.sequence_type == SequenceType.DNA }
    return samples ? samples[0] : null
}

def getDonorDnaSamples(case_record) {
    return case_record.donor_dna_samples.findAll { it.sequence_type == SequenceType.DNA }
}

def getLongitudinalSample(case_record) {
    return case_record.longitudinal_samples ? case_record.longitudinal_samples[0] : null
}

// Sample names
def getTumorDnaSampleName(Map named_args, case_record) {
    def sample
    if (named_args.getOrDefault('primary', false)) {
        sample = getTumorDnaSample(case_record)
    } else {
        sample = getLongitudinalSample(case_record) ?: getTumorDnaSample(case_record)
    }
    return sample?.sample_id
}

def getTumorDnaSampleName(case_record) {
    getTumorDnaSampleName([:], case_record)
}

def getTumorRnaSampleName(case_record) {
    return getTumorRnaSample(case_record)?.sample_id
}

def getNormalDnaSampleName(case_record) {
    return getNormalDnaSample(case_record)?.sample_id
}

def getDonorDnaSampleNames(case_record) {
    return getDonorDnaSamples(case_record).collect { it.sample_id }
}

// Files - Tumor DNA
def getTumorDnaFastq(case_record) {
    return getTumorDnaSample(case_record)?.files?.get(FileType.FASTQ)
}

def getTumorDnaReduxInput(case_record) {
    def d = hasReduxData(getTumorDnaSample(case_record))
    return d
}

def hasTumorDnaFastq(case_record) {
    return getTumorDnaFastq(case_record) != null
}

def hasTumorDnaReduxInput(case_record) {
    return getTumorDnaReduxInput(case_record) != null
}

// Files - Normal DNA
def getNormalDnaFastq(case_record) {
    return getNormalDnaSample(case_record)?.files?.get(FileType.FASTQ)
}

def getNormalDnaReduxInput(case_record) {
    def d = hasReduxData(getNormalDnaSample(case_record))
    return d
}

def hasNormalDnaFastq(case_record) {
    return getNormalDnaFastq(case_record) != null
}

def hasNormalDnaReduxInput(case_record) {
    return getNormalDnaReduxInput(case_record) != null
}

// Files - Donor DNA
def hasDonorDnaFastqs(case_record) {
    return getDonorDnaSamples(case_record).any { it.files?.get(FileType.FASTQ) != null }
}

// Files - Tumor RNA
def getTumorRnaFastq(case_record) {
    return getTumorRnaSample(case_record)?.files?.get(FileType.FASTQ)
}

def hasTumorRnaFastq(case_record) {
    return getTumorRnaFastq(case_record) != null
}

// Status
def hasTumorDna(case_record) {
    return hasInput(getTumorDnaSample(case_record), FileType.ALN) || hasTumorDnaReduxInput(case_record) || hasTumorDnaFastq(case_record)
}

def hasNormalDna(case_record) {
    return hasInput(getNormalDnaSample(case_record), FileType.ALN) || hasNormalDnaReduxInput(case_record) || hasNormalDnaFastq(case_record)
}

def hasNormalDnaAlignment(case_record) {
    return hasInput(getNormalDnaSample(case_record), FileType.ALN) || hasNormalDnaReduxInput(case_record)
}

def hasTumorRna(case_record) {
    return hasInput(getTumorRnaSample(case_record), FileType.ALN) || hasTumorRnaFastq(case_record)
}

def hasReduxData(sample) {
    if (! sample) {
        return null
    }
    return sample.files.get(FileType.ALN_REDUX) ?: sample.files.get(FileType.REDUX_DIR)
}

def getInput(sample, filetype) {
    return sample?.files?.get(filetype)
}

def hasInput(sample, filetype) {
    return getInput(sample, filetype) != null
}

// REDUX alignment and index retrieval
def getTumorReduxDirAlignment(case_record, redux_dir) {
    return getReduxDirAlignment(getTumorDnaSampleName(case_record), redux_dir)
}

def getNormalReduxDirAlignment(case_record, redux_dir) {
    return getReduxDirAlignment(getNormalDnaSampleName(case_record), redux_dir)
}

def getDonorReduxDirAlignments(case_record, donor_dirs) {
    def dirs = donor_dirs == null ? [] : (donor_dirs instanceof List ? donor_dirs : [donor_dirs])
    return getDonorDnaSamples(case_record).collect { donor ->
        def dir = dirs.find { it?.name == "redux_${donor.sample_id}" }
        getReduxDirAlignment(donor.sample_id, dir)
    }
}

def getReduxDirAlignment(sample_name, redux_dir) {
    if (! redux_dir) {
        return [null, null]
    }

    def redux_cram = redux_dir.resolve("${sample_name}.redux.cram")
    if (redux_cram.exists()) {
        return [redux_cram, nextflow.Nextflow.file("${redux_cram.toUriString()}.crai")]
    }

    def redux_bam = redux_dir.resolve("${sample_name}.redux.bam")
    return [redux_bam, nextflow.Nextflow.file("${redux_bam.toUriString()}.bai")]
}

// REDUX TSV retrieval
def getTumorReduxTsvs(case_record, redux_dir) {
    return getReduxTsvs(getTumorDnaSampleName(case_record), redux_dir)
}

def getNormalReduxTsvs(case_record, redux_dir) {
    return getReduxTsvs(getNormalDnaSampleName(case_record), redux_dir)
}

def getDonorReduxTsvs(case_record, donor_dirs) {
    def dirs = donor_dirs == null ? [] : (donor_dirs instanceof List ? donor_dirs : [donor_dirs])
    return getDonorDnaSamples(case_record).collectMany { donor ->
        def dir = dirs.find { it?.name == "redux_${donor.sample_id}" }
        getReduxTsvs(donor.sample_id, dir)
    }
}

def getReduxTsvs(sample_name, redux_dir) {

    if (! redux_dir) {
        return []
    }

    def tsv_exts = [
        'redux.bqr.tsv',
        'redux.jitter_params.tsv',
        'redux.ms_table.tsv.gz',

        'redux.duplicate_freq.tsv',
        'redux.msi_prediction.tsv',

        'umi_coord_freq.tsv',
        'umi_edit_distance.tsv',
        'umi_nucleotide_freq.tsv',
    ]

    return tsv_exts.collect { ext -> redux_dir.resolve("${sample_name}.${ext}") }.findAll { it.exists() }
}

// Existing alignment (plan ALN, else REDUX alignment) to feed REDUX
def getTumorDnaAln(case_record) {
    return getInput(getTumorDnaSample(case_record), FileType.ALN) ?: getInput(getTumorDnaSample(case_record), FileType.ALN_REDUX)
}

def getNormalDnaAln(case_record) {
    return getInput(getNormalDnaSample(case_record), FileType.ALN) ?: getInput(getNormalDnaSample(case_record), FileType.ALN_REDUX)
}

def hasTumorDnaAln(case_record) { return getTumorDnaAln(case_record) != null }
def hasNormalDnaAln(case_record) { return getNormalDnaAln(case_record) != null }
def hasDonorDnaAlns(case_record) {
    return getDonorDnaSamples(case_record).any { getInput(it, FileType.ALN) != null || getInput(it, FileType.ALN_REDUX) != null }
}
