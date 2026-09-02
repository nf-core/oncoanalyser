//
// Sample accessors for the nf-core/oncoanalyser pipeline
//

include { FileType } from './types_enums'

// Generic primitives
def getInput(sample, filetype) {
    return sample?.files?.get(filetype)
}

def hasInput(sample, filetype) {
    return getInput(sample, filetype) != null
}

def hasAnyInput(sample, filetypes) {
    return filetypes.any { hasInput(sample, it) }
}

def firstOrNull(list) {
    return list ? list.first() : null
}

// Sample selectors
def getSamples(case_record) {
    return case_record.normal_dna_samples + case_record.donor_dna_samples + case_record.tumor_dna_samples + case_record.tumor_rna_samples + case_record.longitudinal_samples
}

def getSamples(case_record, sampletypes, sequencetypes) {
    return getSamples(case_record).findAll { s -> sampletypes.contains(s.sample_type) && sequencetypes.contains(s.sequence_type) }
}

def getTumorDnaSample(case_record) {
    return firstOrNull(case_record.tumor_dna_samples)
}

def getNormalDnaSample(case_record) {
    return firstOrNull(case_record.normal_dna_samples)
}

def getDonorDnaSample(case_record) {
    return firstOrNull(case_record.donor_dna_samples)
}

def getDonorDnaSamples(case_record) {
    return case_record.donor_dna_samples
}

def getTumorRnaSample(case_record) {
    return firstOrNull(case_record.tumor_rna_samples)
}

def getLongitudinalSample(case_record) {
    return firstOrNull(case_record.longitudinal_samples)
}

// Sample names
def getTumorDnaSampleName(case_record) {
    return getTumorDnaSample(case_record)?.sample_id
}

def getNormalDnaSampleName(case_record) {
    return getNormalDnaSample(case_record)?.sample_id
}

def getDonorDnaSampleNames(case_record) {
    return getDonorDnaSamples(case_record).collect { it.sample_id }
}

def getTumorRnaSampleName(case_record) {
    return getTumorRnaSample(case_record)?.sample_id
}

def getLongitudinalSampleName(case_record) {
    return getLongitudinalSample(case_record)?.sample_id
}

// Predicates
def hasTumorDna(case_record) {
    return hasAnyInput(getTumorDnaSample(case_record), [FileType.ALN, FileType.ALN_REDUX, FileType.REDUX_DIR, FileType.FASTQ])
}

def hasNormalDna(case_record) {
    return hasAnyInput(getNormalDnaSample(case_record), [FileType.ALN, FileType.ALN_REDUX, FileType.REDUX_DIR, FileType.FASTQ])
}

def hasNormalDnaAlignment(case_record) {
    return hasAnyInput(getNormalDnaSample(case_record), [FileType.ALN, FileType.ALN_REDUX, FileType.REDUX_DIR])
}

def hasTumorRna(case_record) {
    return hasAnyInput(getTumorRnaSample(case_record), [FileType.ALN, FileType.FASTQ])
}

def hasTumorDnaFastq(case_record) {
    return hasInput(getTumorDnaSample(case_record), FileType.FASTQ)
}

def hasNormalDnaFastq(case_record) {
    return hasInput(getNormalDnaSample(case_record), FileType.FASTQ)
}

def hasTumorRnaFastq(case_record) {
    return hasInput(getTumorRnaSample(case_record), FileType.FASTQ)
}

def hasDonorDnaFastqs(case_record) {
    return getDonorDnaSamples(case_record).any { hasInput(it, FileType.FASTQ) }
}

def hasAlignmentInput(sample) {
    return hasAnyInput(sample, [FileType.FASTQ, FileType.ALN, FileType.ALN_REDUX, FileType.REDUX_DIR])
}
