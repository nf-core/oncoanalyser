//
// Named record types for the nf-core/oncoanalyser pipeline case / sample model
//

nextflow.enable.types = true

include { SampleType   } from './types'
include { SequenceType } from './types'

record FastqFile {
    read_fwd: Path
    read_rev: Path?
    single_end: Boolean
    library_id: String
    lane: String
    flowcell: String?
    rg_fields: Map
}

record Sample {
    sample_id: String
    case_id: String
    patient_id: String
    sample_type: SampleType
    sequence_type: SequenceType
    files: Map
    generate_redux_tsvs_only: Boolean
}

record Case {
    case_id: String
    patient_id: String
    cancer_type: String?
    normal_dna_samples: List
    donor_dna_samples: List
    tumor_dna_samples: List
    tumor_rna_samples: List
    longitudinal_samples: List
}
