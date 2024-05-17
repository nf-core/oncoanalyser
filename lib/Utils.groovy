//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

import nextflow.Nextflow
import nextflow.splitter.SplitterEx

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
                        System.exit(1)
                    } else {
                        meta.subject_id = it.subject_id
                    }

                    // Sample type
                    def sample_type_enum = Utils.getEnumFromString(it.sample_type, Constants.SampleType)
                    if (!sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.SampleType).join('\n  - ')
                        log.error "received invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        System.exit(1)
                    }

                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (!sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "received invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        System.exit(1)
                    }

                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (!filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "received invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        System.exit(1)
                    }

                    def sample_key = [sample_type_enum, sequence_type_enum]
                    def meta_sample = meta.get(sample_key, [sample_id: it.sample_id])

                    if (meta_sample.sample_id != it.sample_id) {
                        log.error "got unexpected sample name for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${it.sample_id}"
                        System.exit(1)
                    }

                    if (meta_sample.containsKey(filetype_enum) & filetype_enum != Constants.FileType.FASTQ) {
                        log.error "got duplicate file for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${filetype_enum}"
                        System.exit(1)
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
                                    System.exit(1)
                                }

                                if (info_data.containsKey(info_field_enum)) {
                                    log.error "got duplicate info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${info_field_enum}"
                                    System.exit(1)
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
                            System.exit(1)
                        }

                        if (!info_data.containsKey(Constants.InfoField.LANE)) {
                            log.error "missing 'lane' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            System.exit(1)
                        }

                        def (fwd, rev) = it.filepath.tokenize(';')
                        def fastq_key = [info_data[Constants.InfoField.LIBRARY_ID], info_data[Constants.InfoField.LANE]]

                        if (meta_sample.containsKey(fastq_key)) {
                            log.error "got duplicate lane + library_id data for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${fastq_key}"
                            System.exit(1)
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
                            index_str = 'bai'
                        } else if (key === Constants.FileType.BAM_MARKDUPS) {
                            index_enum = Constants.FileType.BAI
                            index_str = 'bai'
                        } else if (key === Constants.FileType.GRIDSS_VCF) {
                            index_enum = Constants.FileType.GRIDSS_VCF_TBI
                            index_str = 'tbi'
                        } else if (key === Constants.FileType.GRIPSS_VCF) {
                            index_enum = Constants.FileType.GRIPSS_VCF_TBI
                            index_str = 'tbi'
                        } else if (key === Constants.FileType.GRIPSS_UNFILTERED_VCF) {
                            index_enum = Constants.FileType.GRIPSS_UNFILTERED_VCF_TBI
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
                            System.exit(1)
                        }

                        meta[sample_key][index_enum] = index_fp

                    }
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
            params.ref_data_virusbreakenddb_path,
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
                    !meta[key].containsKey(Constants.FileType.BAM_MARKDUPS) &&
                    !meta[key].containsKey(Constants.FileType.FASTQ)) {

                    log.error "no BAMs nor BAM_MARKDUPs nor FASTQs provided for ${meta.group_id} ${sample_type}/${sequence_type}\n\n" +
                        "NB: BAMs or BAM_MARKDUPs or FASTQs are always required as they are the basis to determine input sample type."
                    System.exit(1)
                }

            }

            // Apply some required restrictions to targeted mode
            if (run_config.mode === Constants.RunMode.TARGETED) {

                // Do not allow normal DNA
                if (Utils.hasNormalDna(meta)) {
                    log.error "targeted mode is not compatible with the normal DNA BAM provided for ${meta.group_id}\n\n" +
                        "The targeted workflow supports only tumor DNA BAMs (and tumor RNA BAMs for TSO500)"
                    System.exit(1)
                }

                // Do not allow only tumor RNA
                if (Utils.hasTumorRnaBam(meta) && !Utils.hasTumorDna(meta)) {
                    log.error "targeted mode is not compatible with only tumor RNA provided for ${meta.group_id}\n\n" +
                        "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                        "the configured panel."
                    System.exit(1)
                }

                // Restrict tumor RNA inputs to the TSO500 panel
                if (Utils.hasTumorRnaBam(meta) && run_config.panel != 'tso500') {
                    def panel = run_config.panel.toUpperCase()
                        "Only the TSO500 panel supports tumor RNA analysis"
                    System.exit(1)
                }

            }

            // Do not allow normal DNA only
            if (Utils.hasNormalDna(meta) && !Utils.hasTumorDna(meta)) {
                log.error "germline only mode not supported, found only a normal DNA BAM for ${meta.group_id}\n"
                System.exit(1)
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
                log.error "duplicate sample names found for ${meta.group_id}:\n\n" +
                    "${duplicate_message_strs.join("\n")}"
                System.exit(1)
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
            System.exit(1)
        }

        // Refuse to create STAR index for reference genome containing ALTs, refer to Slack channel
        def run_star_index = run_config.stages.alignment && run_config.has_rna_fastq && !params.ref_data_genome_star_index

        if (run_star_index && has_alt_contigs) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to create the STAR index for a reference genome with ALT contigs.\n" +
                "  Please review https://github.com/alexdobin/STAR docs or contact us on Slack.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        // Require that an input GTF file is provided when creating STAR index
        if (run_star_index && !params.ref_data_genome_gtf) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Creating a STAR index requires the appropriate genome transcript annotations\n" +
                "  as a GTF file. Please contact us on Slack for further information."
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
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
            log.error "recieved an invalid run mode: '${run_mode}'. Valid options are:\n  - ${run_modes_str}"
            System.exit(1)
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


    // Files
    static public getTumorDnaFastq(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    static public getTumorDnaBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    static public getTumorDnaMarkdupsBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAM_MARKDUPS, null)
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

    static public hasTumorDnaMarkdupsBam(meta) {
        return getTumorDnaMarkdupsBam(meta) !== null
    }


    static public getNormalDnaFastq(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    static public getNormalDnaBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    static public getNormalDnaMarkdupsBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAM_MARKDUPS, null)
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

    static public hasNormalDnaMarkdupsBam(meta) {
        return getNormalDnaMarkdupsBam(meta) !== null
    }


    static public hasDnaFastq(meta) {
        return hasNormalDnaFastq(meta) || hasTumorDnaFastq(meta)
    }

    static public hasDnaMarkdupsBam(meta) {
        return hasNormalDnaMarkdupsBam(meta) || hasTumorDnaMarkdupsBam(meta)
    }


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
        return hasTumorDnaBam(meta) || hasTumorDnaMarkdupsBam(meta) || hasTumorDnaFastq(meta)
    }

    static public hasNormalDna(meta) {
        return hasNormalDnaBam(meta) || hasNormalDnaMarkdupsBam(meta) || hasNormalDnaFastq(meta)
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
