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

    def sample_keys = [
        [Constants.SampleType.TUMOR, Constants.SequenceType.DNA],
        [Constants.SampleType.TUMOR, Constants.SequenceType.RNA],
        [Constants.SampleType.NORMAL, Constants.SequenceType.DNA],
    ]

    inputs.each { meta ->

        // Require ALN or ALN_REDUX or REDUX_DIR or FASTQs for each defined sample type
        // NOTE(SW): repeating key pairs above to avoid having to duplicate error messages
        sample_keys.each { key ->

            if (! meta.containsKey(key)) {
                return
            }

            def (sample_type, sequence_type) = key

            if (
                ! meta[key].containsKey(Constants.FileType.FASTQ) &&
                ! meta[key].containsKey(Constants.FileType.ALN) &&
                ! meta[key].containsKey(Constants.FileType.ALN_REDUX) &&
                ! meta[key].containsKey(Constants.FileType.REDUX_DIR)
            ) {
                log.error "no alignments (or REDUX alignments / directory) nor FASTQ files provided for ${meta.group_id} ${sample_type}/${sequence_type}\n\n" +
                    "NB: At least one of these files is required as they are the basis to determine input sample type."
                exit 1
            }

        }

        // Do not allow donor sample without normal sample
        if (hasDonorDna(meta) && ! hasNormalDna(meta)) {
            log.error "a donor sample but no normal sample was found for ${meta.group_id}\n\n" +
                "Analysis with a donor sample requires a normal sample."
            exit 1
        }

        // Apply some required restrictions to targeted mode
        if (run_config.mode == Constants.RunMode.TARGETED) {

            // Do not allow donor DNA
            if (hasDonorDna(meta)) {
                log.error "targeted mode is not compatible with the donor DNA BAM/CRAM provided for ${meta.group_id}\n\n" +
                    "The targeted workflow supports only tumor and normal DNA BAM/CRAMs (and tumor RNA BAM/CRAMs for TSO500)"
                exit 1
            }

            // Do not allow only tumor RNA
            if (hasTumorRna(meta) && ! hasTumorDna(meta)) {
                log.error "targeted mode is not compatible with only tumor RNA provided for ${meta.group_id}\n\n" +
                    "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                    "the configured panel."
                exit 1
            }

        }

        // Do not allow normal DNA only
        if (hasNormalDna(meta) && ! hasTumorDna(meta)) {
            log.error "found only normal DNA input for ${meta.group_id} but germline only analysis is not supported"
            exit 1
        }

        // Enforce unique samples names within groups
        def sample_ids_duplicated = sample_keys
            .groupBy { meta.getOrDefault(it, [:]).getOrDefault('sample_id', null) }
            .findResults { k, v -> k != null && v.size() > 1 ? [k, v] : null }

        if (sample_ids_duplicated) {
            def duplicate_message_strs = sample_ids_duplicated.collect { sample_id, keys ->
                def key_strs = keys.collect { sample_type, sequence_type -> "${sample_type}/${sequence_type}" }
                return "  * ${sample_id}: ${key_strs.join(", ")}"
            }
            log.error "duplicate sample names found for ${meta.group_id}:\n\n${duplicate_message_strs.join("\n")}"
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


// Sample records
def getTumorDnaSample(meta) {
    return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
}

def getTumorRnaSample(meta) {
    return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.RNA], [:])
}

def getNormalDnaSample(meta) {
    return meta.getOrDefault([Constants.SampleType.NORMAL, Constants.SequenceType.DNA], [:])
}

def getDonorDnaSample(meta) {
    return meta.getOrDefault([Constants.SampleType.DONOR, Constants.SequenceType.DNA], [:])
}

// Sample names
def getTumorDnaSampleName(Map named_args, meta) {
    def meta_sample = getTumorDnaSample(meta)
    def sample_id

    if (named_args.getOrDefault('primary', false)) {
        sample_id = meta_sample['sample_id']
    } else {
        sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])
    }

    return sample_id
}

def getTumorDnaSampleName(meta) {
    getTumorDnaSampleName([:], meta)
}

def getTumorRnaSampleName(meta) {
    return getTumorRnaSample(meta)['sample_id']
}

def getNormalDnaSampleName(meta) {
    return getNormalDnaSample(meta)['sample_id']
}

def getDonorDnaSampleName(meta) {
    return getDonorDnaSample(meta)['sample_id']
}


// Files - Tumor DNA
def getTumorDnaFastq(meta) {
    return getTumorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
}

def getTumorDnaBam(meta) {
    return getTumorDnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
}

def getTumorDnaReduxInput(meta) {
    def meta_sample = getTumorDnaSample(meta)
    return hasReduxData(meta_sample) ?: null
}

def getTumorDnaBai(meta) {
    return getTumorDnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
}


def hasTumorDnaFastq(meta) {
    return getTumorDnaFastq(meta) != null
}

def hasTumorDnaBam(meta) {
    return getTumorDnaBam(meta) != null
}

def hasTumorDnaReduxInput(meta) {
    return getTumorDnaReduxInput(meta) != null
}


// Files - Normal DNA
def getNormalDnaFastq(meta) {
    return getNormalDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
}

def getNormalDnaBam(meta) {
    return getNormalDnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
}

def getNormalDnaReduxInput(meta) {
    def meta_sample = getNormalDnaSample(meta)
    return hasReduxData(meta_sample) ?: null
}

def getNormalDnaBai(meta) {
    return getNormalDnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
}


def hasNormalDnaFastq(meta) {
    return getNormalDnaFastq(meta) != null
}

def hasNormalDnaBam(meta) {
    return getNormalDnaBam(meta) != null
}

def hasNormalDnaReduxInput(meta) {
    return getNormalDnaReduxInput(meta) != null
}

def hasDnaFastq(meta) {
    return hasNormalDnaFastq(meta) || hasTumorDnaFastq(meta)
}

def hasDnaReduxInput(meta) {
    return hasNormalDnaReduxInput(meta) || hasTumorDnaReduxInput(meta)
}


// Files - Donor DNA
def getDonorDnaFastq(meta) {
    return getDonorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
}

def getDonorDnaBam(meta) {
    return getDonorDnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
}

def getDonorDnaReduxInput(meta) {
    def meta_sample = getDonorDnaSample(meta)
    return hasReduxData(meta_sample) ?: null
}

def getDonorDnaBai(meta) {
    return getDonorDnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
}


def hasDonorDnaFastq(meta) {
    return getDonorDnaFastq(meta) != null
}

def hasDonorDnaBam(meta) {
    return getDonorDnaBam(meta) != null
}

def hasDonorDnaReduxInput(meta) {
    return getDonorDnaReduxInput(meta) != null
}


// Files - Tumor RNA
def getTumorRnaFastq(meta) {
    return getTumorRnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
}

def getTumorRnaBam(meta) {
    return getTumorRnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
}

def getTumorRnaBai(meta) {
    return getTumorRnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
}


def hasTumorRnaFastq(meta) {
    return getTumorRnaFastq(meta) != null
}

def hasTumorRnaBam(meta) {
    return getTumorRnaBam(meta) != null
}


// Status
def hasTumorDna(meta) {
    return hasTumorDnaBam(meta) || hasTumorDnaReduxInput(meta) || hasTumorDnaFastq(meta)
}

def hasNormalDna(meta) {
    return hasNormalDnaBam(meta) || hasNormalDnaReduxInput(meta) || hasNormalDnaFastq(meta)
}

def hasDonorDna(meta) {
    return hasDonorDnaBam(meta) || hasDonorDnaReduxInput(meta) || hasDonorDnaFastq(meta)
}

def hasTumorRna(meta) {
    return hasTumorRnaBam(meta) || hasTumorRnaFastq(meta)
}

def hasReduxData(meta_sample) {
    return meta_sample.getOrDefault(Constants.FileType.ALN_REDUX, null) || meta_sample.getOrDefault(Constants.FileType.REDUX_DIR, null)
}


// REDUX alignment and index retrieval
def getTumorReduxDirAlignment(meta, redux_dir) {
    return getReduxDirAlignment(getTumorDnaSampleName(meta), redux_dir)
}

def getNormalReduxDirAlignment(meta, redux_dir) {
    return getReduxDirAlignment(getNormalDnaSampleName(meta), redux_dir)
}

def getDonorReduxDirAlignment(meta, redux_dir) {
    return getReduxDirAlignment(getDonorDnaSampleName(meta), redux_dir)
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
def getTumorReduxTsvs(meta, redux_dir) {
    return getReduxTsvs(getTumorDnaSampleName(meta), redux_dir)
}

def getNormalReduxTsvs(meta, redux_dir) {
    return getReduxTsvs(getNormalDnaSampleName(meta), redux_dir)
}

def getDonorReduxTsvs(meta, redux_dir) {
    return getReduxTsvs(getDonorDnaSampleName(meta), redux_dir)
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


// Misc
def getInput(meta, key_set) {
    def keys = key_set.combinations().collect { filetype, sampletype, sequencetype -> return [filetype, [sampletype, sequencetype]] }

    def result = []
    def matched = keys.find { key ->
        def (key_filetype, key_sample) = key
        meta.containsKey(key_sample) && meta[key_sample].containsKey(key_filetype)
    }
    if (matched) {
        def (key_filetype, key_sample) = matched
        result = meta[key_sample].get(key_filetype)
    }
    return result
}

def hasExistingInput(meta, key) {
    return getInput(meta, key) != []
}

def selectCurrentOrExisting(val, meta, key) {
    if (hasExistingInput(meta, key)) {
        return getInput(meta, key)
    } else {
        return val
    }
}

def getDnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ meta ]
    def ch_inputs_tumor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = hasExistingInput(meta, Constants.INPUT.ALN_DNA_TUMOR)
            runnable: hasTumorDnaFastq(meta) && ! has_existing
            skip: true
        }

    def ch_inputs_normal_sorted = ch_inputs
        .branch { meta ->
            def has_existing = hasExistingInput(meta, Constants.INPUT.ALN_DNA_NORMAL)
            runnable: hasNormalDnaFastq(meta) && ! has_existing
            skip: true
        }

    def ch_inputs_donor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = hasExistingInput(meta, Constants.INPUT.ALN_DNA_DONOR)
            runnable: hasDonorDnaFastq(meta) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta -> [meta, getTumorDnaSample(meta), 'tumor'] },
            ch_inputs_normal_sorted.runnable.map { meta -> [meta, getNormalDnaSample(meta), 'normal'] },
            ch_inputs_donor_sorted.runnable.map { meta -> [meta, getDonorDnaSample(meta), 'donor'] },
        )
        .flatMap { meta, meta_sample, sample_type ->
            meta_sample
                .getAt(Constants.FileType.FASTQ)
                .collect { key, d ->
                    def (library_id, lane, flowcell) = key

                    def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

                    def fastq_info = [
                        'sample_id': sample_id,
                        'library_id': library_id,
                        'lane': lane,
                        'sample_type': sample_type,
                        'rg_fields': d.rg_fields,
                    ]

                    if (flowcell) {
                         fastq_info.flowcell = flowcell
                    }

                    return [meta, fastq_info, d['fwd'], d['rev']]
                }
        }

    return channel.empty()
        .mix(
            ch_fastqs,
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, [:], [], []] },
            ch_inputs_normal_sorted.skip.map { meta -> [meta, [:], [], []] },
            ch_inputs_donor_sorted.skip.map { meta -> [meta, [:], [], []] },
        )
}

def getRnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ meta ]
    def ch_inputs_sorted = ch_inputs
        .branch { meta ->
            def has_existing = hasExistingInput(meta, Constants.INPUT.ALN_RNA_TUMOR)
            runnable: hasTumorRnaFastq(meta) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = ch_inputs_sorted.runnable
        .flatMap { meta ->
            def meta_sample = getTumorRnaSample(meta)
            meta_sample
                .getAt(Constants.FileType.FASTQ)
                .collect { key, d ->
                    def (library_id, lane, flowcell) = key

                    def fastq_info = [
                        'sample_id': meta_sample.sample_id,
                        'library_id': library_id,
                        'lane': lane,
                        'rg_fields': d.rg_fields,
                    ]

                    if (flowcell) {
                         fastq_info.flowcell = flowcell
                    }

                    return [meta, fastq_info, d['fwd'], d['rev']]
                }
        }

    return channel.empty()
        .mix(
            ch_fastqs,
            ch_inputs_sorted.skip.map { meta -> [meta, [:], [], []] },
        )
}
