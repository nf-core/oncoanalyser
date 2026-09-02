//
// Alignment and REDUX accessors for the nf-core/oncoanalyser pipeline
//

include { FileType } from './types_enums'
include { getDonorDnaSamples } from './accessors_samples'
include { getInput } from './accessors_samples'
include { getNormalDnaSample } from './accessors_samples'
include { getNormalDnaSampleName } from './accessors_samples'
include { getTumorDnaSample } from './accessors_samples'
include { getTumorDnaSampleName } from './accessors_samples'
include { hasAnyInput } from './accessors_samples'

def getTumorDnaAln(case_record) {
    return getInput(getTumorDnaSample(case_record), FileType.ALN) ?: getInput(getTumorDnaSample(case_record), FileType.ALN_REDUX)
}

def getNormalDnaAln(case_record) {
    return getInput(getNormalDnaSample(case_record), FileType.ALN) ?: getInput(getNormalDnaSample(case_record), FileType.ALN_REDUX)
}

def hasTumorDnaAln(case_record) { return getTumorDnaAln(case_record) != null }
def hasNormalDnaAln(case_record) { return getNormalDnaAln(case_record) != null }
def hasDonorDnaAlns(case_record) {
    return getDonorDnaSamples(case_record).any { hasAnyInput(it, [FileType.ALN, FileType.ALN_REDUX]) }
}

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
    if (! redux_bam.exists()) {
        error "no REDUX alignment (.redux.bam or .redux.cram) found for ${sample_name} in ${redux_dir}"
    }
    return [redux_bam, nextflow.Nextflow.file("${redux_bam.toUriString()}.bai")]
}

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
