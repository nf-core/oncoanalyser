//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import nextflow.Nextflow

class Utils {

    public static parseInput(input_fp_str, stub_run, log) {

        if (! input_fp_str) {
            log.error "Missing required --input argument"
            Nextflow.exit(1)
        }

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
                    if (! sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.SampleType).join('\n  - ')
                        log.error "received invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        Nextflow.exit(1)
                    }

                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (! sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "received invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        Nextflow.exit(1)
                    }

                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (! filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "received invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        Nextflow.exit(1)
                    }

                    def sample_key = [sample_type_enum, sequence_type_enum]
                    def meta_sample = meta.get(sample_key, [:])

                    // Info data
                    def info_data = [:]
                    if (it.containsKey('info')) {
                        // Parse
                        it.info
                            .tokenize(';')
                            .each { e ->

                                def k
                                def v
                                if (e.contains(':')) {
                                   (k, v) = e.split(':', 2)
                                } else {
                                   k = e
                                }

                                def info_field_enum = Utils.getEnumFromString(k, Constants.InfoField)

                                if (! info_field_enum) {
                                    def info_field_str = Utils.getEnumNames(Constants.InfoField).join('\n  - ')
                                    log.error "received invalid info field: '${k}'. Valid options are:\n  - ${info_field_str}"
                                    Nextflow.exit(1)
                                }

                                if (info_data.containsKey(info_field_enum)) {
                                    log.error "got duplicate info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${info_field_enum}"
                                    Nextflow.exit(1)
                                }

                                if (! v && info_field_enum != Constants.InfoField.LONGITUDINAL_SAMPLE && info_field_enum != Constants.InfoField.GENERATE_REDUX_TSVS_ONLY) {
                                    log.error "got empty value for ${group_id} ${sample_type_enum}/${sequence_type_enum} ${info_field_enum}"
                                    Nextflow.exit(1)
                                }

                                info_data[info_field_enum] = v
                            }

                        // Process
                        if (info_data.containsKey(Constants.InfoField.CANCER_TYPE)) {
                            meta[Constants.InfoField.CANCER_TYPE] = info_data[Constants.InfoField.CANCER_TYPE]
                        }

                        if (info_data.containsKey(Constants.InfoField.GENERATE_REDUX_TSVS_ONLY)) {
                            meta_sample[Constants.InfoField.GENERATE_REDUX_TSVS_ONLY] = true
                        }

                        // Only allow READ_GROUP_OVERRIDES for FASTQ
                        if (filetype_enum != Constants.FileType.FASTQ && info_data.containsKey(Constants.InfoField.READ_GROUP_OVERRIDES)) {
                            log.error "The read_group info field is only applicable to FASTQ input but got '${it.filetype}' for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                    }

                    if (info_data.containsKey(Constants.InfoField.LONGITUDINAL_SAMPLE)) {

                        if (meta_sample.containsKey('longitudinal_sample_id') && meta_sample.longitudinal_sample_id != it.sample_id) {
                            log.error "got multiple longitudinal samples for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${it.sample_id}"
                            Nextflow.exit(1)
                        }

                        meta_sample.longitudinal_sample_id = it.sample_id

                    } else if (meta_sample.containsKey('sample_id') && meta_sample.sample_id != it.sample_id) {

                        log.error "got unexpected sample name for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${it.sample_id}"
                        Nextflow.exit(1)

                    } else {

                        meta_sample.sample_id = it.sample_id

                    }

                    // Disallow AMBER, COBALT, SAGE_APPEND inputs for longitudinal samples; these would clash with primary inputs
                    if (info_data.containsKey(Constants.InfoField.LONGITUDINAL_SAMPLE)) {
                        def longitudinal_disallowed_input_list = [Constants.FileType.AMBER_DIR, Constants.FileType.COBALT_DIR, Constants.FileType.SAGE_APPEND_DIR]
                        def disallowed_inputs = longitudinal_disallowed_input_list.findAll { e -> e == filetype_enum }
                        if (disallowed_inputs) {
                            log.error "got disallowed ${filetype_enum} input for longitudinal sample ${group_id} ${meta_sample.sample_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }
                    }

                    // Filetype uniqueness
                    if (meta_sample.containsKey(filetype_enum) && filetype_enum != Constants.FileType.FASTQ) {
                        log.error "got duplicate file for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${filetype_enum}"
                        Nextflow.exit(1)
                    }

                    // Handle inputs appropriately
                    if (filetype_enum == Constants.FileType.FASTQ) {

                        if (! info_data.containsKey(Constants.InfoField.LIBRARY_ID)) {
                            log.error "missing 'library_id' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                        if (! info_data.containsKey(Constants.InfoField.LANE)) {
                            log.error "missing 'lane' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                        def fastq_entries = it.filepath.tokenize(';')

                        if (fastq_entries.size() != 2) {
                            log.error "expected exactly 2 FASTQ files delimited by ';' (i.e. '<fwd>;<rev>') but found ${fastq_entries.size} " +
                                " files for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                        def (fwd, rev) = fastq_entries
                        def fastq_key = [info_data[Constants.InfoField.LIBRARY_ID], info_data[Constants.InfoField.LANE], info_data.getOrDefault(Constants.InfoField.FLOWCELL, null)]

                        if (! meta_sample.containsKey(filetype_enum)) {
                            meta_sample[filetype_enum] = [:]
                        }

                        if (meta_sample[filetype_enum].containsKey(fastq_key)) {
                            log.error "got duplicate lane + library_id data for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${fastq_key}"
                            Nextflow.exit(1)
                        }

                        def rg_fields = [:]
                        if (info_data.containsKey(Constants.InfoField.READ_GROUP_OVERRIDES)) {
                            rg_fields = Utils.parse_read_group_info(info_data[Constants.InfoField.READ_GROUP_OVERRIDES], log)
                        }

                        meta_sample[filetype_enum][fastq_key] = ['fwd': Utils.getFileObject(fwd), 'rev': Utils.getFileObject(rev), 'rg_fields': rg_fields]

                    } else {

                        meta_sample[filetype_enum] = Utils.getFileObject(it.filepath)

                    }

                    // Record sample key to simplify iteration later on
                    sample_keys << sample_key
                }

                // Set BAM / CRAM as generic alignments, along with their index
                sample_keys.each { sample_key ->

                    def meta_sample = meta[sample_key]

                    if (meta_sample.containsKey(Constants.FileType.BAM)) {
                        meta_sample[Constants.FileType.ALN] = meta_sample.remove(Constants.FileType.BAM)
                    }

                    if (meta_sample.containsKey(Constants.FileType.CRAM)) {
                        meta_sample[Constants.FileType.ALN] = meta_sample.remove(Constants.FileType.CRAM)
                    }

                    if (meta_sample.containsKey(Constants.FileType.BAM_REDUX)) {
                        meta_sample[Constants.FileType.ALN_REDUX] = meta_sample.remove(Constants.FileType.BAM_REDUX)
                    }

                    if (meta_sample.containsKey(Constants.FileType.CRAM_REDUX)) {
                        meta_sample[Constants.FileType.ALN_REDUX] = meta_sample.remove(Constants.FileType.CRAM_REDUX)
                    }

                    if (meta_sample.containsKey(Constants.FileType.BAI)) {
                        meta_sample[Constants.FileType.IDX] = meta_sample.remove(Constants.FileType.BAI)
                    }

                    if (meta_sample.containsKey(Constants.FileType.CRAI)) {
                        meta_sample[Constants.FileType.IDX] = meta_sample.remove(Constants.FileType.CRAI)
                    }

                }

                // Handle REDUX alignments; require TSVs to be present, and promote to directory when files for this sample are found.
                // Allow multiple samples to be in the same dir (e.g. both tumor and normal sample).
                // Existing helpers (getReduxDirAlignment, getReduxTsvs) already resolve files by exact {sample_id}.redux.* filename.
                sample_keys.each { sample_key ->

                    def meta_sample = meta[sample_key]
                    def is_primary = ! meta_sample.containsKey('longitudinal_sample_id')
                    def sample_id = is_primary ? meta_sample.sample_id : meta_sample.longitudinal_sample_id

                    if (stub_run) {
                        return
                    }

                    if (meta_sample.containsKey(Constants.FileType.ALN_REDUX) && meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {
                        log.error "expected either REDUX directory or REDUX alignment but got both for ${meta.group_id} ${sample_id}"
                        Nextflow.exit(1)
                    }

                    def redux_input
                    def redux_aln
                    def redux_dir
                    if (meta_sample.containsKey(Constants.FileType.ALN_REDUX)) {

                        redux_input = redux_aln = meta_sample[Constants.FileType.ALN_REDUX]
                        redux_dir = redux_aln.parent
                        if (! redux_input.isFile()) {
                            log.error "didn't receive file as REDUX alignment for ${meta.group_id} ${sample_id}: ${redux_input}"
                            Nextflow.exit(1)
                        }

                    } else if (meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {

                        redux_input = redux_dir = meta_sample[Constants.FileType.REDUX_DIR]
                        redux_aln = getReduxDirAlignment(sample_id, redux_dir)[0]
                        if (! redux_input.isDirectory()) {
                            log.error "didn't receive directory as REDUX directory for ${meta.group_id} ${sample_id}: ${redux_input}"
                            Nextflow.exit(1)
                        }

                    } else {

                        return

                    }

                    def has_bqr_tsv = redux_dir.resolve("${sample_id}.redux.bqr.tsv").exists()
                    def has_jitter_tsv = redux_dir.resolve("${sample_id}.redux.jitter_params.tsv").exists()
                    def has_ms_tsv = redux_dir.resolve("${sample_id}.redux.ms_table.tsv.gz").exists()
                    def has_redux_tsvs = has_bqr_tsv && has_jitter_tsv && has_ms_tsv

                    def has_colocated_index = nextflow.Nextflow.file("${redux_aln.toUriString()}.bai").exists() || nextflow.Nextflow.file("${redux_aln.toUriString()}.crai").exists()

                    def generate_tsvs_only = meta_sample.getOrDefault(Constants.InfoField.GENERATE_REDUX_TSVS_ONLY, false)

                    if (! has_redux_tsvs && ! generate_tsvs_only) {

                        log.error "no REDUX TSVs provided or found for ${meta.group_id} ${sample_id}: ${redux_input}"
                        Nextflow.exit(1)

                    }

                    if (! has_colocated_index) {

                        if (meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {
                            log.error "required index not located in REDUX directory: ${meta.group_id} ${sample_id}: ${redux_input}"
                            Nextflow.exit(1)
                        }

                    }

                    if (has_colocated_index && ! meta_sample.containsKey(Constants.FileType.REDUX_DIR)) {
                        meta_sample.remove(Constants.FileType.ALN_REDUX)
                        meta_sample[Constants.FileType.REDUX_DIR] = redux_dir
                    }

                    if (! has_redux_tsvs && generate_tsvs_only) {
                        meta_sample.remove(Constants.FileType.REDUX_DIR)
                        meta_sample[Constants.FileType.ALN_REDUX] = redux_aln
                    }

                    if (meta_sample.containsKey(Constants.FileType.ALN_REDUX) && meta_sample.containsKey(Constants.FileType.IDX) && ! has_redux_tsvs && ! generate_tsvs_only) {
                        log.error "REDUX alignments without colocated TSVs requires generate_redux_tsvs_only to be set in the samplesheet: ${meta.group_id} ${sample_id}: ${redux_input}"
                        Nextflow.exit(1)
                    }

                    if (meta_sample.containsKey(Constants.FileType.REDUX_DIR) && has_redux_tsvs && generate_tsvs_only) {
                        log.warn "REDUX directory already contains TSVs, ignoring generate_redux_tsvs_only flag for: ${meta.group_id} ${sample_id}: ${redux_input}"
                    }

                }

                // Check that required indexes are provided or are accessible
                sample_keys.each { sample_key ->

                    meta[sample_key]*.key.each { key ->

                        def meta_sample = meta[sample_key]
                        def is_primary = ! meta_sample.containsKey('longitudinal_sample_id')
                        def sample_id = is_primary ? meta_sample.sample_id : meta_sample.longitudinal_sample_id

                        if (stub_run) {
                            return
                        }

                        def aln
                        if (key == Constants.FileType.ALN || key == Constants.FileType.ALN_REDUX) {
                            aln = meta_sample[key]
                        } else if (key == Constants.FileType.REDUX_DIR) {
                            def d = getReduxDirAlignment(sample_id, meta_sample[key])
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
                            def (sample_type, sequence_type) = sample_key
                            log.error "got bad alignment type for ${meta.group_id} ${sample_type}/${sequence_type}: ${key}: ${aln}"
                            Nextflow.exit(1)
                        }

                        if (meta_sample.containsKey(Constants.FileType.IDX) && meta_sample[Constants.FileType.IDX].name.endsWith(index_ext)) {
                            return
                        }

                        def fp = aln.toUriString()
                        def index_fp = nextflow.Nextflow.file("${fp}.${index_ext}")

                        if (! index_fp.exists() && ! stub_run) {
                            def (sample_type, sequence_type) = sample_key
                            log.error "no index provided or found for ${meta.group_id} ${sample_type}/${sequence_type}: ${key}: ${fp}"
                            Nextflow.exit(1)
                        }

                        if (key == Constants.FileType.ALN || key == Constants.FileType.ALN_REDUX) {
                            meta_sample[Constants.FileType.IDX] = index_fp
                        }
                    }
                }

                // For purity estimation with WISP, require primary normal DNA BAM when an AMBER directory is provided
                def meta_tumor_dna = meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
                def longitudinal = meta_tumor_dna.containsKey('longitudinal_sample_id')
                def has_amber_dir = meta_tumor_dna.containsKey(Constants.FileType.AMBER_DIR)
                def has_normal_dna_bam = Utils.hasNormalDnaBam(meta) || Utils.hasNormalDnaReduxInput(meta)

                if (longitudinal && has_amber_dir && ! has_normal_dna_bam) {
                    log.error "AMBER input was provided without the required primary normal DNA BAM for ${meta.group_id}"
                    Nextflow.exit(1)
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

            def fp = Utils.getFileObject(fp_str)

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

    public static void validateInput(inputs, run_config, params, log) {

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
                    Nextflow.exit(1)
                }

            }

            // Do not allow donor sample without normal sample
            if (Utils.hasDonorDna(meta) && ! Utils.hasNormalDna(meta)) {
                log.error "a donor sample but no normal sample was found for ${meta.group_id}\n\n" +
                    "Analysis with a donor sample requires a normal sample."
                Nextflow.exit(1)
            }

            // Apply some required restrictions to targeted mode
            if (run_config.mode == Constants.RunMode.TARGETED) {

                // Do not allow donor DNA
                if (Utils.hasDonorDna(meta)) {
                    log.error "targeted mode is not compatible with the donor DNA BAM/CRAM provided for ${meta.group_id}\n\n" +
                        "The targeted workflow supports only tumor and normal DNA BAM/CRAMs (and tumor RNA BAM/CRAMs for TSO500)"
                    Nextflow.exit(1)
                }

                // Do not allow only tumor RNA
                if (Utils.hasTumorRna(meta) && ! Utils.hasTumorDna(meta)) {
                    log.error "targeted mode is not compatible with only tumor RNA provided for ${meta.group_id}\n\n" +
                        "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                        "the configured panel."
                    Nextflow.exit(1)
                }

            }

            // Do not allow normal DNA only
            if (Utils.hasNormalDna(meta) && ! Utils.hasTumorDna(meta)) {
                log.error "found only normal DNA input for ${meta.group_id} but germline only analysis is not supported"
                Nextflow.exit(1)
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
                Nextflow.exit(1)
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
            Nextflow.exit(1)
        }

        // Refuse to create STAR index for reference genome containing ALTs, refer to Slack channel
        def run_star_index = run_config.stages.alignment && run_config.has_rna_fastq && ! params.ref_data_genome_star_index

        if (run_star_index && has_alt_contigs) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to create the STAR index for a reference genome with ALT contigs.\n" +
                "  Please review https://github.com/alexdobin/STAR docs or contact us on Slack.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Require that an input GTF file is provided when creating STAR index
        if (run_star_index && ! params.ref_data_genome_gtf) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Creating a STAR index requires the appropriate genome transcript annotations\n" +
                "  as a GTF file. Please contact us on Slack for further information.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Require --isofox_gene_ids argument to be provided in PANEL_RESOURCE_CREATION when RNA inputs are present
        if (run_config.mode == Constants.RunMode.PANEL_RESOURCE_CREATION && run_config.has_rna && ! params.isofox_gene_ids) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Running the panel resource creation workflow with RNA requires that the\n" +
                "  --isofox_gene_ids argument is set with an appropriate input file.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

    }

    public static parse_read_group_info(rg_info_raw, log) {
        def escape_char = "\u0000"
        def validate_rg_tags = ['BC', 'CN', 'DS', 'DT', 'FO', 'ID', 'KS', 'LB', 'PG', 'PI', 'PL', 'PM', 'PU', 'SM']

        def fields = [:]
        def rg_info_escaped = rg_info_raw.replace('||', escape_char)
        rg_info_escaped.split('\\|').each { field_str_escaped ->
            def field_str = field_str_escaped.replace(escape_char, '|')
            if (! field_str.contains('=')) {
                log.error "Received bad read group field (must be in format `<name>=<value>`): ${field_str}"
                Nextflow.exit(1)
            }

            def (name, value) = field_str.split('=', 2)
            if (! validate_rg_tags.contains(name)) {
                log.error "Received bad read group tag '${name}' in: ${rg_info_raw}"
                Nextflow.exit(1)
            }

            if (! value) {
                log.error "Received empty read group value for '${name}' in: ${rg_info_raw}"
                Nextflow.exit(1)
            }

            fields[name] = value
        }

        return fields
    }

    public static getSequencingPlatformPons(hmf_data, sequencing_platform_string, log) {
        def sequencing_platform = Utils.getEnumFromString(sequencing_platform_string, Constants.SequencingPlatform)
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
                Nextflow.exit(1)
            }
        }
    }

    public static getEnumFromString(s, e) {
        try {
            return e.valueOf(s.toUpperCase())
        } catch (java.lang.IllegalArgumentException err) {
            return null
        }
    }

    public static getEnumNames(e) {
        e
            .values()
            *.name()
            *.toLowerCase()
    }


    public static getFileObject(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    public static getRunMode(run_mode, log) {
        def run_mode_enum = Utils.getEnumFromString(run_mode, Constants.RunMode)
        if (! run_mode_enum) {
            def run_modes_str = Utils.getEnumNames(Constants.RunMode).join('\n  - ')
            log.error "received an invalid run mode: '${run_mode}'. Valid options are:\n  - ${run_modes_str}"
            Nextflow.exit(1)
        }
        return run_mode_enum
    }


    // Sample records
    public static getTumorDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
    }

    public static getTumorRnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.RNA], [:])
    }

    public static getNormalDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.NORMAL, Constants.SequenceType.DNA], [:])
    }

    public static getDonorDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.DONOR, Constants.SequenceType.DNA], [:])
    }

    // Sample names
    public static getTumorDnaSampleName(Map named_args, meta) {
        def meta_sample = getTumorDnaSample(meta)
        def sample_id

        if (named_args.getOrDefault('primary', false)) {
            sample_id = meta_sample['sample_id']
        } else {
            sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])
        }

        return sample_id
    }

    public static getTumorDnaSampleName(meta) {
        getTumorDnaSampleName([:], meta)
    }

    public static getTumorRnaSampleName(meta) {
        return getTumorRnaSample(meta)['sample_id']
    }

    public static getNormalDnaSampleName(meta) {
        return getNormalDnaSample(meta)['sample_id']
    }

    public static getDonorDnaSampleName(meta) {
        return getDonorDnaSample(meta)['sample_id']
    }


    // Files - Tumor DNA
    public static getTumorDnaFastq(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getTumorDnaBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
    }

    public static getTumorDnaReduxInput(meta) {
        def meta_sample = getTumorDnaSample(meta)
        return hasReduxData(meta_sample) ?: null
    }

    public static getTumorDnaBai(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
    }


    public static hasTumorDnaFastq(meta) {
        return getTumorDnaFastq(meta) != null
    }

    public static hasTumorDnaBam(meta) {
        return getTumorDnaBam(meta) != null
    }

    public static hasTumorDnaReduxInput(meta) {
        return getTumorDnaReduxInput(meta) != null
    }


    // Files - Normal DNA
    public static getNormalDnaFastq(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getNormalDnaBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
    }

    public static getNormalDnaReduxInput(meta) {
        def meta_sample = getNormalDnaSample(meta)
        return hasReduxData(meta_sample) ?: null
    }

    public static getNormalDnaBai(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
    }


    public static hasNormalDnaFastq(meta) {
        return getNormalDnaFastq(meta) != null
    }

    public static hasNormalDnaBam(meta) {
        return getNormalDnaBam(meta) != null
    }

    public static hasNormalDnaReduxInput(meta) {
        return getNormalDnaReduxInput(meta) != null
    }

    public static hasDnaFastq(meta) {
        return hasNormalDnaFastq(meta) || hasTumorDnaFastq(meta)
    }

    public static hasDnaReduxInput(meta) {
        return hasNormalDnaReduxInput(meta) || hasTumorDnaReduxInput(meta)
    }


    // Files - Donor DNA
    public static getDonorDnaFastq(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getDonorDnaBam(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
    }

    public static getDonorDnaReduxInput(meta) {
        def meta_sample = getDonorDnaSample(meta)
        return hasReduxData(meta_sample) ?: null
    }

    public static getDonorDnaBai(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
    }


    public static hasDonorDnaFastq(meta) {
        return getDonorDnaFastq(meta) != null
    }

    public static hasDonorDnaBam(meta) {
        return getDonorDnaBam(meta) != null
    }

    public static hasDonorDnaReduxInput(meta) {
        return getDonorDnaReduxInput(meta) != null
    }


    // Files - Tumor RNA
    public static getTumorRnaFastq(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getTumorRnaBam(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.ALN, null)
    }

    public static getTumorRnaBai(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.IDX, null)
    }


    public static hasTumorRnaFastq(meta) {
        return getTumorRnaFastq(meta) != null
    }

    public static hasTumorRnaBam(meta) {
        return getTumorRnaBam(meta) != null
    }


    // Status
    public static hasTumorDna(meta) {
        return hasTumorDnaBam(meta) || hasTumorDnaReduxInput(meta) || hasTumorDnaFastq(meta)
    }

    public static hasNormalDna(meta) {
        return hasNormalDnaBam(meta) || hasNormalDnaReduxInput(meta) || hasNormalDnaFastq(meta)
    }

    public static hasDonorDna(meta) {
        return hasDonorDnaBam(meta) || hasDonorDnaReduxInput(meta) || hasDonorDnaFastq(meta)
    }

    public static hasTumorRna(meta) {
        return hasTumorRnaBam(meta) || hasTumorRnaFastq(meta)
    }

    public static hasReduxData(meta_sample) {
        return meta_sample.getOrDefault(Constants.FileType.ALN_REDUX, null) || meta_sample.getOrDefault(Constants.FileType.REDUX_DIR, null)
    }


    // REDUX alignment and index retrieval
    public static getTumorReduxDirAlignment(meta, redux_dir) {
        return getReduxDirAlignment(getTumorDnaSampleName(meta), redux_dir)
    }

    public static getNormalReduxDirAlignment(meta, redux_dir) {
        return getReduxDirAlignment(getNormalDnaSampleName(meta), redux_dir)
    }

    public static getDonorReduxDirAlignment(meta, redux_dir) {
        return getReduxDirAlignment(getDonorDnaSampleName(meta), redux_dir)
    }

    public static getReduxDirAlignment(sample_name, redux_dir) {
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
    public static getTumorReduxTsvs(meta, redux_dir) {
        return getReduxTsvs(Utils.getTumorDnaSampleName(meta), redux_dir)
    }

    public static getNormalReduxTsvs(meta, redux_dir) {
        return getReduxTsvs(Utils.getNormalDnaSampleName(meta), redux_dir)
    }

    public static getDonorReduxTsvs(meta, redux_dir) {
        return getReduxTsvs(Utils.getDonorDnaSampleName(meta), redux_dir)
    }

    public static getReduxTsvs(sample_name, redux_dir) {

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
    public static getInput(meta, key_set) {
        def keys = key_set
            .combinations()
            .collect { filetype, sampletype, sequencetype -> return [filetype, [sampletype, sequencetype]] }

        def result = []
        for (key in keys) {
            def (key_filetype, key_sample) = key
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
