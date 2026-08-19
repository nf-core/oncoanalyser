//
// Sample parsing, validation, and accessor helpers for the nf-core/oncoanalyser pipeline
//

def createStubPlaceholders(params) {

    def fps = [
        params.ref_data_genome_alt,
        params.ref_data_genome_bwamem2_index,
        params.ref_data_genome_dict,
        params.ref_data_genome_fai,
        params.ref_data_genome_fasta,
        params.ref_data_genome_gridss_index,
        params.ref_data_genome_gtf,
        params.ref_data_genome_star_index,
    ]

    params.hmf_data_paths[params.genome_version.toString()]
        .each { k, v ->
            fps << "${params.ref_data_hmf_data_path.replaceAll('/$', '')}/${v}"
        }

    if (params.panel != null) {
        params.panel_data_paths[params.panel.toLowerCase()][params.genome_version.toString()]
            .each { k, v ->
                fps << "${params.ref_data_panel_data_path.replaceAll('/$', '')}/${v}"
            }
    }

    fps.each { fp_str ->
        if (fp_str == null) {
            return
        }

        def fp = getFileObject(fp_str)

        if (! fp_str || fp.exists()) {
            return
        }

        if (fp_str.endsWith('/')) {
            fp.mkdirs()
        } else {
            fp.getParent().mkdirs()
            fp.toFile().createNewFile()
        }
    }

}

def validateInput(inputs, run_config, params, log) {

    inputs.each { case_record ->

        // Require ALN or ALN_REDUX or REDUX_DIR or FASTQs for each defined sample
        getSamples(case_record).each { sample ->
            if (! hasAlignmentInput(sample)) {
                log.error "no alignments (or REDUX alignments / directory) nor FASTQ files provided for ${case_record.case_id} ${sample.sample_type}/${sample.sequence_type}\n\n" +
                    "NB: At least one of these files is required as they are the basis to determine input sample type."
                exit 1
            }
        }

        // Do not allow donor sample without normal sample
        if (case_record.donor_dna_samples && ! case_record.normal_dna_samples) {
            log.error "a donor sample but no normal sample was found for ${case_record.case_id}\n\n" +
                "Analysis with a donor sample requires a normal sample."
            exit 1
        }

        // Longitudinal samples require a primary normal DNA alignment when AMBER input is provided
        if (case_record.longitudinal_samples && case_record.directories.containsKey(Constants.FileType.AMBER_DIR) && ! hasNormalDnaAlignment(case_record)) {
            log.error "AMBER input was provided without the required primary normal DNA BAM for ${case_record.case_id}"
            exit 1
        }

        // Apply some required restrictions to targeted mode
        if (run_config.mode == Constants.RunMode.TARGETED) {

            // Do not allow donor DNA
            if (case_record.donor_dna_samples) {
                log.error "targeted mode is not compatible with the donor DNA BAM/CRAM provided for ${case_record.case_id}\n\n" +
                    "The targeted workflow supports only tumor and normal DNA BAM/CRAMs (and tumor RNA BAM/CRAMs for TSO500)"
                exit 1
            }

            // Do not allow only tumor RNA
            if (hasTumorRna(case_record) && ! hasTumorDna(case_record)) {
                log.error "targeted mode is not compatible with only tumor RNA provided for ${case_record.case_id}\n\n" +
                    "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                    "the configured panel."
                exit 1
            }

        }

        // Do not allow normal DNA only
        if (hasNormalDna(case_record) && ! hasTumorDna(case_record)) {
            log.error "found only normal DNA input for ${case_record.case_id} but germline only analysis is not supported"
            exit 1
        }

        // Enforce unique samples names within cases
        def sample_ids_duplicated = getSamples(case_record)
            .groupBy { it.sample_id }
            .findAll { k, v -> v.size() > 1 }

        if (sample_ids_duplicated) {
            def duplicate_message_strs = sample_ids_duplicated.collect { sample_id, samples ->
                def key_strs = samples.collect { "${it.sample_type}/${it.sequence_type}" }
                return "  * ${sample_id}: ${key_strs.join(", ")}"
            }
            log.error "duplicate sample names found for ${case_record.case_id}:\n\n${duplicate_message_strs.join("\n")}"
            exit 1
        }

    }


    // NOTE(SW): the following final config checks are performed here since they require additional information
    // regarding processes that are run and also inputs

    def has_alt_contigs = params.genome_type == 'alt'

    // Ensure that custom genomes with ALT contigs that need indexes built have the required .alt file
    def has_bwa_indexes = (params.ref_data_genome_bwamem2_index && params.ref_data_genome_gridss_index)
    def has_alt_file = params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt
    def run_bwa_or_gridss_index = run_config.stages.alignment && run_config.has_dna_fastq && ! has_bwa_indexes

    if (run_bwa_or_gridss_index && has_alt_contigs && ! has_alt_file) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  The genome .alt file is required when building bwa-mem2 or GRIDSS indexes\n" +
            "  for reference genomes containing ALT contigs\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    // Refuse to create STAR index for reference genome containing ALTs, refer to Slack channel
    def run_star_index = run_config.stages.alignment && run_config.has_rna_fastq && ! params.ref_data_genome_star_index

    if (run_star_index && has_alt_contigs) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Refusing to create the STAR index for a reference genome with ALT contigs.\n" +
            "  Please review https://github.com/alexdobin/STAR docs or contact us on Slack.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    // Require that an input GTF file is provided when creating STAR index
    if (run_star_index && ! params.ref_data_genome_gtf) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Creating a STAR index requires the appropriate genome transcript annotations\n" +
            "  as a GTF file. Please contact us on Slack for further information.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    // Require --isofox_gene_ids argument to be provided in PANEL_RESOURCE_CREATION when RNA inputs are present
    if (run_config.mode == Constants.RunMode.PANEL_RESOURCE_CREATION && run_config.has_rna && ! params.isofox_gene_ids) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Running the panel resource creation workflow with RNA requires that the\n" +
            "  --isofox_gene_ids argument is set with an appropriate input file.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
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

def getSequencingPlatformPons(hmf_data, sequencing_platform_string, log) {
    def sequencing_platform = getEnumFromString(sequencing_platform_string, Constants.SequencingPlatform)
    hmf_data.map { d ->
        if (sequencing_platform == Constants.SequencingPlatform.ILLUMINA) {
            return [
                'esvee_breakends': d.esvee_pon_breakends_illumina,
                'esvee_breakpoints': d.esvee_pon_breakpoints_illumina,
                'sage': d.sage_pon_illumina,
            ]
        } else if (sequencing_platform == Constants.SequencingPlatform.SBX) {
            return [
                'esvee_breakends': d.esvee_pon_breakends_sbx,
                'esvee_breakpoints': d.esvee_pon_breakpoints_sbx,
                'sage': d.sage_pon_sbx,
            ]
        } else if (sequencing_platform == Constants.SequencingPlatform.ULTIMA) {
            return [
                'esvee_breakends': d.esvee_pon_breakends_ultima,
                'esvee_breakpoints': d.esvee_pon_breakpoints_ultima,
                'sage': d.sage_pon_ultima,
            ]
        } else {
            log.error "Got bad sequencing platform: ${sequencing_platform}"
            exit 1
        }
    }
}

def getEnumFromString(s, e) {
    try {
        return e.valueOf(s.toUpperCase())
    } catch (java.lang.IllegalArgumentException err) {
        return null
    }
}

def getEnumNames(e) {
    e
        .values()
        *.name()
        *.toLowerCase()
}

def getEnumFromStringOrFail(value, enum_class, label, log) {
    def enum_value = getEnumFromString(value, enum_class)
    if (! enum_value) {
        def options = getEnumNames(enum_class).join('\n  - ')
        log.error "received invalid ${label}: '${value}'. Valid options are:\n  - ${options}"
        exit 1
    }
    return enum_value
}


def getFileObject(path) {
    return path ? nextflow.Nextflow.file(path) : []
}

def getRunMode(run_mode, log) {
    def run_mode_enum = getEnumFromString(run_mode, Constants.RunMode)
    if (! run_mode_enum) {
        def run_modes_str = getEnumNames(Constants.RunMode).join('\n  - ')
        log.error "received an invalid run mode: '${run_mode}'. Valid options are:\n  - ${run_modes_str}"
        exit 1
    }
    return run_mode_enum
}


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
    return files.containsKey(Constants.FileType.FASTQ) ||
        files.containsKey(Constants.FileType.ALN) ||
        files.containsKey(Constants.FileType.ALN_REDUX) ||
        files.containsKey(Constants.FileType.REDUX_DIR)
}


// Sample records (singular: first matching sample, or null)
def getTumorDnaSample(case_record) {
    def samples = case_record.tumor_dna_samples.findAll { it.sequence_type == Constants.SequenceType.DNA }
    return samples ? samples[0] : null
}

def getTumorRnaSample(case_record) {
    return case_record.tumor_rna_samples ? case_record.tumor_rna_samples[0] : null
}

def getNormalDnaSample(case_record) {
    def samples = case_record.normal_dna_samples.findAll { it.sequence_type == Constants.SequenceType.DNA }
    return samples ? samples[0] : null
}

def getDonorDnaSample(case_record) {
    def samples = case_record.donor_dna_samples.findAll { it.sequence_type == Constants.SequenceType.DNA }
    return samples ? samples[0] : null
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

def getDonorDnaSampleName(case_record) {
    return getDonorDnaSample(case_record)?.sample_id
}

def getLongitudinalSampleName(case_record) {
    return getLongitudinalSample(case_record)?.sample_id
}


// Files - Tumor DNA
def getTumorDnaFastq(case_record) {
    return getTumorDnaSample(case_record)?.files?.get(Constants.FileType.FASTQ)
}

def getTumorDnaBam(case_record) {
    return getTumorDnaSample(case_record)?.files?.get(Constants.FileType.ALN)?.path
}

def getTumorDnaReduxInput(case_record) {
    def d = hasReduxData(getTumorDnaSample(case_record))
    return d ? d.path : null
}

def getTumorDnaBai(case_record) {
    return getTumorDnaSample(case_record)?.files?.get(Constants.FileType.IDX)?.path
}


def hasTumorDnaFastq(case_record) {
    return getTumorDnaFastq(case_record) != null
}

def hasTumorDnaBam(case_record) {
    return getTumorDnaBam(case_record) != null
}

def hasTumorDnaReduxInput(case_record) {
    return getTumorDnaReduxInput(case_record) != null
}


// Files - Normal DNA
def getNormalDnaFastq(case_record) {
    return getNormalDnaSample(case_record)?.files?.get(Constants.FileType.FASTQ)
}

def getNormalDnaBam(case_record) {
    return getNormalDnaSample(case_record)?.files?.get(Constants.FileType.ALN)?.path
}

def getNormalDnaReduxInput(case_record) {
    def d = hasReduxData(getNormalDnaSample(case_record))
    return d ? d.path : null
}

def getNormalDnaBai(case_record) {
    return getNormalDnaSample(case_record)?.files?.get(Constants.FileType.IDX)?.path
}


def hasNormalDnaFastq(case_record) {
    return getNormalDnaFastq(case_record) != null
}

def hasNormalDnaBam(case_record) {
    return getNormalDnaBam(case_record) != null
}

def hasNormalDnaReduxInput(case_record) {
    return getNormalDnaReduxInput(case_record) != null
}

def hasDnaFastq(case_record) {
    return hasNormalDnaFastq(case_record) || hasTumorDnaFastq(case_record)
}

def hasDnaReduxInput(case_record) {
    return hasNormalDnaReduxInput(case_record) || hasTumorDnaReduxInput(case_record)
}


// Files - Donor DNA
def getDonorDnaFastq(case_record) {
    return getDonorDnaSample(case_record)?.files?.get(Constants.FileType.FASTQ)
}

def getDonorDnaBam(case_record) {
    return getDonorDnaSample(case_record)?.files?.get(Constants.FileType.ALN)?.path
}

def getDonorDnaReduxInput(case_record) {
    def d = hasReduxData(getDonorDnaSample(case_record))
    return d ? d.path : null
}

def getDonorDnaBai(case_record) {
    return getDonorDnaSample(case_record)?.files?.get(Constants.FileType.IDX)?.path
}


def hasDonorDnaFastq(case_record) {
    return getDonorDnaFastq(case_record) != null
}

def hasDonorDnaBam(case_record) {
    return getDonorDnaBam(case_record) != null
}

def hasDonorDnaReduxInput(case_record) {
    return getDonorDnaReduxInput(case_record) != null
}


// Files - Tumor RNA
def getTumorRnaFastq(case_record) {
    return getTumorRnaSample(case_record)?.files?.get(Constants.FileType.FASTQ)
}

def getTumorRnaBam(case_record) {
    return getTumorRnaSample(case_record)?.files?.get(Constants.FileType.ALN)?.path
}

def getTumorRnaBai(case_record) {
    return getTumorRnaSample(case_record)?.files?.get(Constants.FileType.IDX)?.path
}


def hasTumorRnaFastq(case_record) {
    return getTumorRnaFastq(case_record) != null
}

def hasTumorRnaBam(case_record) {
    return getTumorRnaBam(case_record) != null
}


// Status
def hasTumorDna(case_record) {
    return hasTumorDnaBam(case_record) || hasTumorDnaReduxInput(case_record) || hasTumorDnaFastq(case_record)
}

def hasNormalDna(case_record) {
    return hasNormalDnaBam(case_record) || hasNormalDnaReduxInput(case_record) || hasNormalDnaFastq(case_record)
}

def hasNormalDnaAlignment(case_record) {
    return hasNormalDnaBam(case_record) || hasNormalDnaReduxInput(case_record)
}

def hasDonorDna(case_record) {
    return hasDonorDnaBam(case_record) || hasDonorDnaReduxInput(case_record) || hasDonorDnaFastq(case_record)
}

def hasTumorRna(case_record) {
    return hasTumorRnaBam(case_record) || hasTumorRnaFastq(case_record)
}

def hasReduxData(sample) {
    if (! sample) {
        return null
    }
    return sample.files.get(Constants.FileType.ALN_REDUX, null) ?: sample.files.get(Constants.FileType.REDUX_DIR, null)
}


// REDUX alignment and index retrieval
def getTumorReduxDirAlignment(case_record, redux_dir) {
    return getReduxDirAlignment(getTumorDnaSampleName(case_record), redux_dir)
}

def getNormalReduxDirAlignment(case_record, redux_dir) {
    return getReduxDirAlignment(getNormalDnaSampleName(case_record), redux_dir)
}

def getDonorReduxDirAlignment(case_record, redux_dir) {
    return getReduxDirAlignment(getDonorDnaSampleName(case_record), redux_dir)
}

def getReduxDirAlignment(sample_name, redux_dir) {
    if (! redux_dir) {
        return [[], []]
    }

    def redux_cram = redux_dir.resolve("${sample_name}.redux.cram")
    if (redux_cram.exists()) {
        return [redux_cram, "${redux_cram.toUriString()}.crai"]
    }

    def redux_bam = redux_dir.resolve("${sample_name}.redux.bam")
    return [redux_bam, "${redux_bam.toUriString()}.bai"]
}


// REDUX TSV retrieval
def getTumorReduxTsvs(case_record, redux_dir) {
    return getReduxTsvs(getTumorDnaSampleName(case_record), redux_dir)
}

def getNormalReduxTsvs(case_record, redux_dir) {
    return getReduxTsvs(getNormalDnaSampleName(case_record), redux_dir)
}

def getDonorReduxTsvs(case_record, redux_dir) {
    return getReduxTsvs(getDonorDnaSampleName(case_record), redux_dir)
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


// Misc: resolve an existing input (or existing output dir) from a case
def getInput(case_record, key_set) {
    def filetypes = key_set[0] instanceof List ? key_set[0] : [key_set[0]]
    def sampletypes = key_set[1] instanceof List ? key_set[1] : [key_set[1]]
    def sequencetypes = key_set[2] instanceof List ? key_set[2] : [key_set[2]]

    // Case-level directories
    def dir_ft = filetypes.find { ft -> case_record.directories.containsKey(ft) }
    if (dir_ft) {
        return case_record.directories[dir_ft].path
    }

    // Sample-level files
    def samples = getSamples(case_record, sampletypes, sequencetypes)
    def match = samples.findResult { sample ->
        def ft = filetypes.find { f -> sample.files.containsKey(f) }
        return ft ? sample.files[ft].path : null
    }
    return match ?: []
}

def hasExistingInput(case_record, key) {
    return getInput(case_record, key) != []
}

def selectCurrentOrExisting(val, case_record, key) {
    if (hasExistingInput(case_record, key)) {
        return getInput(case_record, key)
    } else {
        return val
    }
}

def getDnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ case_record ]
    def ch_inputs_tumor_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasExistingInput(case_record, Constants.INPUT.ALN_DNA_TUMOR)
            runnable: hasTumorDnaFastq(case_record) && ! has_existing
            skip: true
        }

    def ch_inputs_normal_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasExistingInput(case_record, Constants.INPUT.ALN_DNA_NORMAL)
            runnable: hasNormalDnaFastq(case_record) && ! has_existing
            skip: true
        }

    def ch_inputs_donor_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasExistingInput(case_record, Constants.INPUT.ALN_DNA_DONOR)
            runnable: hasDonorDnaFastq(case_record) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ case_record, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { case_record -> [case_record, getTumorDnaSample(case_record), 'tumor'] },
            ch_inputs_normal_sorted.runnable.map { case_record -> [case_record, getNormalDnaSample(case_record), 'normal'] },
            ch_inputs_donor_sorted.runnable.map { case_record -> [case_record, getDonorDnaSample(case_record), 'donor'] },
        )
        .flatMap { case_record, sample, sample_type ->
            sample.files.getAt(Constants.FileType.FASTQ)
                .collect { fastq ->
                    def fastq_info = [
                        'sample_id': sample.sample_id,
                        'library_id': fastq.library_id,
                        'lane': fastq.lane,
                        'sample_type': sample_type,
                        'single_end': fastq.single_end,
                        'rg_fields': fastq.rg_fields,
                    ]

                    if (fastq.flowcell) {
                         fastq_info.flowcell = fastq.flowcell
                    }

                    return [case_record, fastq_info, fastq.read_fwd, fastq.read_rev ?: []]
                }
        }

    return channel.empty()
        .mix(
            ch_fastqs,
            ch_inputs_tumor_sorted.skip.map { case_record -> [case_record, [:], [], []] },
            ch_inputs_normal_sorted.skip.map { case_record -> [case_record, [:], [], []] },
            ch_inputs_donor_sorted.skip.map { case_record -> [case_record, [:], [], []] },
        )
}

def getRnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ case_record ]
    def ch_inputs_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasExistingInput(case_record, Constants.INPUT.ALN_RNA_TUMOR)
            runnable: hasTumorRnaFastq(case_record) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ case_record, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = ch_inputs_sorted.runnable
        .flatMap { case_record ->
            def sample = getTumorRnaSample(case_record)
            sample.files
                .getAt(Constants.FileType.FASTQ)
                .collect { fastq ->
                    def fastq_info = [
                        'sample_id': sample.sample_id,
                        'library_id': fastq.library_id,
                        'lane': fastq.lane,
                        'rg_fields': fastq.rg_fields,
                    ]

                    if (fastq.flowcell) {
                         fastq_info.flowcell = fastq.flowcell
                    }

                    return [case_record, fastq_info, fastq.read_fwd, fastq.read_rev]
                }
        }

    return channel.empty()
        .mix(
            ch_fastqs,
            ch_inputs_sorted.skip.map { case_record -> [case_record, [:], [], []] },
        )
}
