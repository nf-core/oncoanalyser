//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import nextflow.Nextflow

class Utils {

    public static parseInput(input_fp_str, stub_run, log) {

        // NOTE(SW): using NF .splitCsv channel operator, hence should be easily interchangable with NF syntax

        def input_fp = Utils.getFileObject(input_fp_str)
        def inputs = nextflow.splitter.SplitterEx.splitCsv(input_fp, [header: true])
            .groupBy { it['group_id'] }
            .collect { group_id, entries ->

                def meta = [group_id: group_id]
                def sample_keys = [] as Set

                // Process each entry
                entries.each {
                    // Add subject id if absent or check if current matches existing
                    if (meta.containsKey('subject_id') && meta.subject_id != it.subject_id) {
                        log.error "got unexpected subject name for ${group_id} ${meta.subject_id}: ${it.subject_id}"
                        Nextflow.exit(1)
                    } else {
                        meta.subject_id = it.subject_id
                    }

                    // Sample type
                    def sample_type_enum = Utils.getEnumFromString(it.sample_type, Constants.SampleType)
                    if (!sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.SampleType).join('\n  - ')
                        log.error "received invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        Nextflow.exit(1)
                    }

                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (!sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "received invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        Nextflow.exit(1)
                    }

                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (!filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "received invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        Nextflow.exit(1)
                    }

                    def sample_key = [sample_type_enum, sequence_type_enum]
                    def meta_sample = meta.get(sample_key, [sample_id: it.sample_id])

                    if (meta_sample.sample_id != it.sample_id) {
                        log.error "got unexpected sample name for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${it.sample_id}"
                        Nextflow.exit(1)
                    }

                    if (meta_sample.containsKey(filetype_enum) & filetype_enum != Constants.FileType.FASTQ) {
                        log.error "got duplicate file for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${filetype_enum}"
                        Nextflow.exit(1)
                    }

                    // Info data
                    def info_data = [:]
                    if (it.containsKey('info')) {
                        // Parse
                        it.info
                            .tokenize(';')
                            .each { e ->
                                def (k, v) = e.tokenize(':')
                                def info_field_enum = Utils.getEnumFromString(k, Constants.InfoField)

                                if (!info_field_enum) {
                                    def info_field_str = Utils.getEnumNames(Constants.InfoField).join('\n  - ')
                                    log.error "received invalid info field: '${k}'. Valid options are:\n  - ${info_field_str}"
                                    Nextflow.exit(1)
                                }

                                if (info_data.containsKey(info_field_enum)) {
                                    log.error "got duplicate info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${info_field_enum}"
                                    Nextflow.exit(1)
                                }

                                info_data[info_field_enum] = v
                            }

                        // Process
                        if (info_data.containsKey(Constants.InfoField.CANCER_TYPE)) {
                            meta[Constants.InfoField.CANCER_TYPE] = info_data[Constants.InfoField.CANCER_TYPE]
                        }

                    }


                    // Handle inputs appropriately
                    if (filetype_enum === Constants.FileType.FASTQ) {

                        if (!info_data.containsKey(Constants.InfoField.LIBRARY_ID)) {
                            log.error "missing 'library_id' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                        if (!info_data.containsKey(Constants.InfoField.LANE)) {
                            log.error "missing 'lane' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                        def (fwd, rev) = it.filepath.tokenize(';')
                        def fastq_key = [info_data[Constants.InfoField.LIBRARY_ID], info_data[Constants.InfoField.LANE]]

                        if (meta_sample.containsKey(fastq_key)) {
                            log.error "got duplicate lane + library_id data for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${fastq_key}"
                            Nextflow.exit(1)
                        }

                        if (!meta_sample.containsKey(filetype_enum)) {
                            meta_sample[filetype_enum] = [:]
                        }

                        meta_sample[filetype_enum][fastq_key] = ['fwd': fwd, 'rev': rev]

                    } else {

                        meta_sample[filetype_enum] = Utils.getFileObject(it.filepath)

                    }

                    // Record sample key to simplify iteration later on
                    sample_keys << sample_key
                }

                // Check that required indexes are provided or are accessible
                sample_keys.each { sample_key ->

                    meta[sample_key]*.key.each { key ->

                        // NOTE(SW): I was going to use two maps but was unable to get an enum map to compile

                        def index_enum
                        def index_str

                        if (key === Constants.FileType.BAM) {
                            index_enum = Constants.FileType.BAI
                            index_str = (meta[sample_key][key].toString().endsWith('cram')) ? 'crai' : 'bai'
                        } else if (key === Constants.FileType.BAM_REDUX) {
                            index_enum = Constants.FileType.BAI
                            index_str = 'bai'
                        } else if (key === Constants.FileType.ESVEE_VCF) {
                            index_enum = Constants.FileType.ESVEE_VCF_TBI
                            index_str = 'tbi'
                        } else if (key === Constants.FileType.SAGE_VCF) {
                            index_enum = Constants.FileType.SAGE_VCF_TBI
                            index_str = 'tbi'
                        } else {
                            return
                        }

                        if (meta[sample_key].containsKey(index_enum)) {
                            return
                        }

                        def fp = meta[sample_key][key].toUriString()
                        def index_fp = nextflow.Nextflow.file("${fp}.${index_str}")

                        if (!index_fp.exists() && !stub_run) {
                            def (sample_type, sequence_type) = sample_key
                            log.error "no index provided or found for ${meta.group_id} ${sample_type}/${sequence_type}: ${key}: ${fp}"
                            Nextflow.exit(1)
                        }

                        meta[sample_key][index_enum] = index_fp

                    }
                }

                // Check that REDUX TSVs are present
                sample_keys.each { sample_key ->

                    if(stub_run)
                        return

                    def meta_sample = meta[sample_key]
                    def sample_id = meta_sample.sample_id

                    if(!meta_sample.containsKey(Constants.FileType.BAM_REDUX))
                        return

                    if(meta_sample.containsKey(Constants.FileType.BAM)) {
                        log.error "${Constants.FileType.BAM} and ${Constants.FileType.BAM_REDUX} provided for sample ${sample_id}. Please only provide one or the other"
                        Nextflow.exit(1)
                    }

                    def bam_path = meta_sample[Constants.FileType.BAM_REDUX]
                    def bam_dir = bam_path.getParent().toUriString()

                    // Get user specified TSV paths
                    def jitter_tsv   = meta_sample[Constants.FileType.REDUX_JITTER_TSV]
                    def ms_tsv       = meta_sample[Constants.FileType.REDUX_MS_TSV]

                    // If TSV paths not provided, default to TSV paths in the same dir as the BAM
                    jitter_tsv   = jitter_tsv   ?: "${bam_dir}/${sample_id}.jitter_params.tsv"
                    ms_tsv       = ms_tsv       ?: "${bam_dir}/${sample_id}.ms_table.tsv.gz"

                    jitter_tsv   = nextflow.Nextflow.file(jitter_tsv)
                    ms_tsv       = nextflow.Nextflow.file(ms_tsv)

                    def missing_tsvs = [:]
                    if(!jitter_tsv.exists()) missing_tsvs[Constants.FileType.REDUX_JITTER_TSV] = jitter_tsv
                    if(!ms_tsv.exists())     missing_tsvs[Constants.FileType.REDUX_MS_TSV] = ms_tsv

                    if(missing_tsvs.size() > 0){

                        def error_message = []

                        error_message.add("When only specifying filetype ${Constants.FileType.BAM_REDUX} in the sample sheet, make sure the REDUX BAM and TSVs are in the same dir:")
                        error_message.add("${bam_path.toUriString()} (${Constants.FileType.BAM_REDUX})")
                        missing_tsvs.each { error_message.add("${it.value} (missing expected ${it.key})") }
                        error_message.add("")
                        error_message.add("Alternatively, provide the TSV paths in the sample sheet using filetype values: " +
                            "${Constants.FileType.REDUX_JITTER_TSV}, " +
                            "${Constants.FileType.REDUX_MS_TSV}"
                        )

                        log.error error_message.join("\n")
                        Nextflow.exit(1)
                    }

                    // Set parsed REDUX TSV paths in metadata object
                    meta_sample[Constants.FileType.REDUX_JITTER_TSV] = jitter_tsv
                    meta_sample[Constants.FileType.REDUX_MS_TSV] = ms_tsv
                }

                return meta
            }

        return inputs
    }

    public static void createStubPlaceholders(params) {

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

        if(params.panel !== null) {
            params.panel_data_paths[params.panel][params.genome_version.toString()]
                .each { k, v ->
                    fps << "${params.ref_data_panel_data_path.replaceAll('/$', '')}/${v}"
                }
        }

        fps.each { fp_str ->
            if (fp_str === null) return

            def fp = Utils.getFileObject(fp_str)

            if (!fp_str || fp.exists()) return

            if (fp_str.endsWith('/')) {
                fp.mkdirs()
            } else {
                fp.getParent().mkdirs()
                fp.toFile().createNewFile()
            }
        }

    }

    public static void validateInput(inputs, run_config, params, log) {

        def sample_keys = [
            [Constants.SampleType.TUMOR, Constants.SequenceType.DNA],
            [Constants.SampleType.TUMOR, Constants.SequenceType.RNA],
            [Constants.SampleType.NORMAL, Constants.SequenceType.DNA],
        ]

        inputs.each { meta ->

            // Require BAMs or BAM_MARKDUPs or FASTQs for each defined sample type
            // NOTE(SW): repeating key pairs above to avoid having to duplicate error messages
            sample_keys.each { key ->

                if (!meta.containsKey(key)) {
                    return
                }

                def (sample_type, sequence_type) = key

                if (!meta[key].containsKey(Constants.FileType.BAM) &&
                    !meta[key].containsKey(Constants.FileType.BAM_REDUX) &&
                    !meta[key].containsKey(Constants.FileType.FASTQ)) {

                    log.error "no BAMs nor BAM_MARKDUPs nor FASTQs provided for ${meta.group_id} ${sample_type}/${sequence_type}\n\n" +
                        "NB: BAMs or BAM_MARKDUPs or FASTQs are always required as they are the basis to determine input sample type."
                    Nextflow.exit(1)
                }

            }

            // Do not allow donor sample without normal sample
            if (Utils.hasDonorDna(meta) && ! Utils.hasNormalDna(meta)) {
                log.error "a donor sample but not normal sample was found for ${meta.group_id}\n\n" +
                    "Analysis with a donor sample requires a normal sample."
                Nextflow.exit(1)
            }

            // Apply some required restrictions to targeted mode
            if (run_config.mode === Constants.RunMode.TARGETED) {

                // Do not allow donor DNA
                if (Utils.hasDonorDna(meta)) {
                    log.error "targeted mode is not compatible with the donor DNA BAM provided for ${meta.group_id}\n\n" +
                        "The targeted workflow supports only tumor and normal DNA BAMs (and tumor RNA BAMs for TSO500)"
                    Nextflow.exit(1)
                }

                // Do not allow only tumor RNA
                if (Utils.hasTumorRna(meta) && !Utils.hasTumorDna(meta)) {
                    log.error "targeted mode is not compatible with only tumor RNA provided for ${meta.group_id}\n\n" +
                        "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                        "the configured panel."
                    Nextflow.exit(1)
                }

            }

            // Do not allow normal DNA only
            if (Utils.hasNormalDna(meta) && !Utils.hasTumorDna(meta)) {
                log.error "found only normal DNA input for ${meta.group_id} but germline only analysis is not supported"
                Nextflow.exit(1)
            }

            // Do not allow CRAM RNA input
            if (Utils.hasTumorRnaBam(meta) && Utils.getTumorRnaBam(meta).toString().endsWith('cram')) {
                log.error "found tumor RNA CRAM input for ${meta.group_id} but RNA CRAM input is not supported"
                Nextflow.exit(1)
            }

            // Enforce unique samples names within groups
            def sample_ids_duplicated = sample_keys
                .groupBy { meta.getOrDefault(it, [:]).getOrDefault('sample_id', null) }
                .findResults { k, v -> k !== null & v.size() > 1 ? [k, v] : null }

            if (sample_ids_duplicated) {
                def duplicate_message_strs = sample_ids_duplicated.collect { sample_id, keys ->
                    def key_strs = keys.collect { sample_type, sequence_type -> "${sample_type}/${sequence_type}" }
                    return "  * ${sample_id}: ${key_strs.join(", ")}"
                }
                log.error "duplicate sample names found for ${meta.group_id}:\n\n${duplicate_message_strs.join("\n")}"
                Nextflow.exit(1)
            }

        }


        // NOTE(SW): the follwing final config checks are performed here since they require additional information
        // regarding processes that are run and also inputs

        def has_alt_contigs = params.genome_type == 'alt'

        // Ensure that custom genomes with ALT contigs that need indexes built have the required .alt file
        def has_bwa_indexes = (params.ref_data_genome_bwamem2_index && params.ref_data_genome_gridss_index)
        def has_alt_file = params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt
        def run_bwa_or_gridss_index = run_config.stages.alignment && run_config.has_dna_fastq && !has_bwa_indexes

        if (run_bwa_or_gridss_index && has_alt_contigs && !has_alt_file) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  The genome .alt file is required when building bwa-mem2 or GRIDSS indexes\n" +
                "  for reference genomes containing ALT contigs\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Refuse to create STAR index for reference genome containing ALTs, refer to Slack channel
        def run_star_index = run_config.stages.alignment && run_config.has_rna_fastq && !params.ref_data_genome_star_index

        if (run_star_index && has_alt_contigs) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to create the STAR index for a reference genome with ALT contigs.\n" +
                "  Please review https://github.com/alexdobin/STAR docs or contact us on Slack.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Require that an input GTF file is provided when creating STAR index
        if (run_star_index && !params.ref_data_genome_gtf) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Creating a STAR index requires the appropriate genome transcript annotations\n" +
                "  as a GTF file. Please contact us on Slack for further information.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

    }

    static public getEnumFromString(s, e) {
        try {
            return e.valueOf(s.toUpperCase())
        } catch(java.lang.IllegalArgumentException err) {
            return null
        }
    }

    public static getEnumNames(e) {
        e
            .values()
            *.name()
            *.toLowerCase()
    }


    static public getFileObject(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    static public getRunMode(run_mode, log) {
        def run_mode_enum = Utils.getEnumFromString(run_mode, Constants.RunMode)
        if (!run_mode_enum) {
            def run_modes_str = Utils.getEnumNames(Constants.RunMode).join('\n  - ')
            log.error "received an invalid run mode: '${run_mode}'. Valid options are:\n  - ${run_modes_str}"
            Nextflow.exit(1)
        }
        return run_mode_enum
    }


    // Sample records
    static public getTumorDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
    }

    static public getTumorRnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.RNA], [:])
    }

    static public getNormalDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.NORMAL, Constants.SequenceType.DNA], [:])
    }

    static public getDonorDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.DONOR, Constants.SequenceType.DNA], [:])
    }

    // Sample names
    static public getTumorDnaSampleName(meta) {
        return getTumorDnaSample(meta)['sample_id']
    }

    static public getTumorRnaSampleName(meta) {
        return getTumorRnaSample(meta)['sample_id']
    }

    static public getNormalDnaSampleName(meta) {
        return getNormalDnaSample(meta)['sample_id']
    }

    static public getDonorDnaSampleName(meta) {
        return getDonorDnaSample(meta)['sample_id']
    }


    // Files - Tumor DNA
    static public getTumorDnaFastq(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    static public getTumorDnaBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    static public getTumorDnaReduxBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAM_REDUX, null)
    }

    static public getTumorDnaBai(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
    }


    static public hasTumorDnaFastq(meta) {
        return getTumorDnaFastq(meta) !== null
    }

    static public hasTumorDnaBam(meta) {
        return getTumorDnaBam(meta) !== null
    }

    static public hasTumorDnaReduxBam(meta) {
        return getTumorDnaReduxBam(meta) !== null
    }


    // Files - Normal DNA
    static public getNormalDnaFastq(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    static public getNormalDnaBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    static public getNormalDnaReduxBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAM_REDUX, null)
    }
    static public getNormalDnaBai(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
    }


    static public hasNormalDnaFastq(meta) {
        return getNormalDnaFastq(meta) !== null
    }

    static public hasNormalDnaBam(meta) {
        return getNormalDnaBam(meta) !== null
    }

    static public hasNormalDnaReduxBam(meta) {
        return getNormalDnaReduxBam(meta) !== null
    }

    static public hasDnaFastq(meta) {
        return hasNormalDnaFastq(meta) || hasTumorDnaFastq(meta)
    }

    static public hasDnaReduxBam(meta) {
        return hasNormalDnaReduxBam(meta) || hasTumorDnaReduxBam(meta)
    }


    // Files - Donor DNA
    static public getDonorDnaFastq(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    static public getDonorDnaBam(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    static public getDonorDnaReduxBam(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.BAM_REDUX, null)
    }

    static public getDonorDnaBai(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
    }


    static public hasDonorDnaFastq(meta) {
        return getDonorDnaFastq(meta) !== null
    }

    static public hasDonorDnaBam(meta) {
        return getDonorDnaBam(meta) !== null
    }

    static public hasDonorDnaReduxBam(meta) {
        return getDonorDnaReduxBam(meta) !== null
    }


    // Files - Tumor RNA
    static public getTumorRnaFastq(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    static public getTumorRnaBam(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    static public getTumorRnaBai(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
    }


    static public hasTumorRnaFastq(meta) {
        return getTumorRnaFastq(meta) !== null
    }

    static public hasTumorRnaBam(meta) {
        return getTumorRnaBam(meta) !== null
    }


    // Status
    static public hasTumorDna(meta) {
        return hasTumorDnaBam(meta) || hasTumorDnaReduxBam(meta) || hasTumorDnaFastq(meta)
    }

    static public hasNormalDna(meta) {
        return hasNormalDnaBam(meta) || hasNormalDnaReduxBam(meta) || hasNormalDnaFastq(meta)
    }

    static public hasDonorDna(meta) {
        return hasDonorDnaBam(meta) || hasDonorDnaReduxBam(meta) || hasDonorDnaFastq(meta)
    }

    static public hasTumorRna(meta) {
        return hasTumorRnaBam(meta) || hasTumorRnaFastq(meta)
    }


    // Misc
    public static getInput(meta, key) {

        def result = []
        def (key_filetype, key_filetypes, key_sequencetypes) = key

        for (key_sample in [key_filetypes, key_sequencetypes].combinations()) {
            if (meta.containsKey(key_sample) && meta[key_sample].containsKey(key_filetype)) {
                result = meta[key_sample].get(key_filetype)
                break
            }
        }
        return result
    }

    public static hasExistingInput(meta, key) {
        return getInput(meta, key) != []
    }

    public static selectCurrentOrExisting(val, meta, key) {
        if (hasExistingInput(meta, key)) {
            return getInput(meta, key)
        } else {
            return val
        }
    }

}
