//
// Record constructors for the nf-core/oncoanalyser case/sample model
//
// Built with the legacy `record()` factory, which returns a RecordMap (a LinkedHashMap
// subclass carrying a Record marker). Fields are read with dot access or bracket access,
// e.g. `case.case_id` or `case['case_id']`.
//

def FastqFile(read_fwd, read_rev, single_end, library_id, lane, flowcell, rg_fields) {
    return record(
        read_fwd: read_fwd,
        read_rev: read_rev,
        single_end: single_end,
        library_id: library_id,
        lane: lane,
        flowcell: flowcell,
        rg_fields: rg_fields,
    )
}

def SampleRecord(sample_id, case_id, patient_id, sample_type, sequence_type, files, generate_redux_tsvs_only) {
    return record(
        sample_id: sample_id,
        case_id: case_id,
        patient_id: patient_id,
        sample_type: sample_type,
        sequence_type: sequence_type,
        files: files,
        generate_redux_tsvs_only: generate_redux_tsvs_only,
    )
}

def CaseRecord(case_id, patient_id, cancer_type, normal_dna_samples, donor_dna_samples, tumor_dna_samples, tumor_rna_samples, longitudinal_samples, directories) {
    return record(
        case_id: case_id,
        patient_id: patient_id,
        cancer_type: cancer_type,
        normal_dna_samples: normal_dna_samples,
        donor_dna_samples: donor_dna_samples,
        tumor_dna_samples: tumor_dna_samples,
        tumor_rna_samples: tumor_rna_samples,
        longitudinal_samples: longitudinal_samples,
        directories: directories,
    )
}
