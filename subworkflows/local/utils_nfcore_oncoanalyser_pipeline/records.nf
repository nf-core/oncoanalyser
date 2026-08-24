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

//
// Per-process meta records. `key` is the channel-plumbing field used by
// groupByMeta / restoreMeta / joinMeta (accessed via getAt('key')); the rest are
// the fields each process group actually reads off the old `meta` Map.
//

record SampleMeta {
    key: String
    id: String
    sample_id: String
    sample_type: String?
}

record ReferenceMeta {
    key: String
    id: String
    tumor_id: String
    normal_id: String?
    donor_id: String?
}

record NormalReferenceMeta {
    key: String
    id: String
    tumor_id: String
    normal_id: String?
}

record TumorNormalMeta {
    key: String
    id: String
    tumor_id: String?
    normal_id: String?
}

record RnaSampleMeta {
    key: String
    id: String
    sample_id: String
    sample_rna_id: String?
}

record NeoScorerMeta {
    key: String
    id: String
    sample_id: String
    cancer_type: String?
    sample_rna_id: String?
}

record OrangeMeta {
    key: String
    id: String
    tumor_id: String
    cancer_type: String?
    normal_dna_id: String?
    tumor_rna_id: String?
}

record ReduxMeta {
    key: String
    id: String
    sample_id: String
    sample_type: String
}

record WispMeta {
    key: String
    id: String
    patient_id: String
    longitudinal_id: String
    primary_id: String
}

record StarMeta {
    key: String
    id: String
    sample_id: String
}
