//
// Prepare reference data as required
//

include { BWAMEM2_INDEX         } from '../../../modules/nf-core/bwamem2/index/main'
include { BWA_INDEX             } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_DICT         } from '../../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_BWA_INDEX_IMAGE } from '../../../modules/local/gatk4/bwaindeximage/main'
include { STAR_GENOMEGENERATE   } from '../../../modules/nf-core/star/genomegenerate/main'
include { GRIDSS_INDEX          } from '../../../modules/local/gridss/index/main'

include { CUSTOM_EXTRACTTARBALL as DECOMP_BWAMEM2_INDEX } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_GRIDSS_INDEX  } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_HMF_DATA      } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_PANEL_DATA    } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_STAR_INDEX    } from '../../../modules/local/custom/extract_tarball/main'

include { WRITE_REFERENCE_DATA as WRITE_FASTA           } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_FAI             } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_DICT            } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_IMG             } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_BWA_INDEX       } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_GRIDSS_INDEX    } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_STAR_INDEX      } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_HMF_DATA        } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_PANEL_DATA      } from '../../../modules/local/custom/write_reference_data/main'

workflow PREPARE_REFERENCE {
    take:
    prepare_reference_only // boolean: [mandatory] prepare reference only, do not run pipeline
    inputs                 // map:     [optional]  sample metadata
    stages                 // map:     [optional]  processes to run

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    //
    // Determine which resources need to be prepared
    //
    def pipeline_mode = RunModes.Pipeline.fromString(params.mode)

    def prep_config = prepare_reference_only
        ? getConfigForPrepRefOnly(params, log)
        : getConfigForPipelineRun(inputs, pipeline_mode, stages)

    def has_alt_contigs = params.genome_type == refgenome.RefGenomeType.ALT
    def has_alt_file = params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt

    //
    // Set .fasta and main genome indexes, create if required
    //
    ch_genome_version = Channel.value(params.genome_version)

    ch_genome_fasta = Channel.empty()
    if (prep_config.require_fasta) {
        ch_genome_fasta = Channel.fromPath(params.ref_data_genome_fasta)
    }

    ch_genome_fai = Channel.empty()
    if (prep_config.require_fai) {

        ch_genome_fai = getRefFileChannel('ref_data_genome_fai')
        if (!params.ref_data_genome_fai) {
            SAMTOOLS_FAIDX(ch_genome_fasta)
            ch_genome_fai = SAMTOOLS_FAIDX.out.fai
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
    }

    ch_genome_dict = Channel.empty()
    if (prep_config.require_dict) {

        ch_genome_dict = getRefFileChannel('ref_data_genome_dict')
        if (!params.ref_data_genome_dict) {
            SAMTOOLS_DICT(ch_genome_fasta)
            ch_genome_dict = SAMTOOLS_DICT.out.dict
            ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)
        }
    }

    ch_genome_img = Channel.empty()
    if (prep_config.require_img) {

        ch_genome_img = getRefFileChannel('ref_data_genome_img')
        if (!params.ref_data_genome_img) {
            GATK4_BWA_INDEX_IMAGE(ch_genome_fasta)
            ch_genome_img = GATK4_BWA_INDEX_IMAGE.out.img
            ch_versions = ch_versions.mix(GATK4_BWA_INDEX_IMAGE.out.versions)
        }
    }

    //
    // Set bwa-mem2 index, unpack or create if required
    //
    ch_genome_bwamem2_index = Channel.empty()
    if (prep_config.require_bwamem2_index) {

        if (!params.ref_data_genome_bwamem2_index) {

            if(has_alt_contigs && !has_alt_file)
                error "For ref genomes with ALT contigs, an .alt file is required when building bwa-mem2 indexes"

            BWAMEM2_INDEX(
                ch_genome_fasta,
                params.ref_data_genome_alt ? file(params.ref_data_genome_alt) : [],
            )
            ch_genome_bwamem2_index = BWAMEM2_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

        } else if (params.ref_data_genome_bwamem2_index.endsWith('.tar.gz')) {

            ch_genome_bwamem2_index_inputs = Channel.fromPath(params.ref_data_genome_bwamem2_index)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_BWAMEM2_INDEX(ch_genome_bwamem2_index_inputs)
            ch_genome_bwamem2_index = DECOMP_BWAMEM2_INDEX.out.extracted_dir

        } else {

            ch_genome_bwamem2_index = getRefFileChannel('ref_data_genome_bwamem2_index')

        }
    }

    //
    // Set GRIDSS index, unpack or create if required
    //
    ch_genome_gridss_index = Channel.empty()
    if (prep_config.require_gridss_index) {

        if (!params.ref_data_genome_gridss_index) {

            if(has_alt_contigs && !has_alt_file)
                error "For ref genomes with ALT contigs, an .alt file is required when building GRIDSS indexes"

            BWA_INDEX(
                ch_genome_fasta,
                params.ref_data_genome_alt ? file(params.ref_data_genome_alt) : [],
            )
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

            GRIDSS_INDEX(
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_dict,
                BWA_INDEX.out.index,
            )
            ch_genome_gridss_index = GRIDSS_INDEX.out.index
            ch_versions = ch_versions.mix(GRIDSS_INDEX.out.versions)

        } else if (params.ref_data_genome_gridss_index.endsWith('.tar.gz')) {

            ch_genome_gridss_index_inputs = Channel.fromPath(params.ref_data_genome_gridss_index)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_GRIDSS_INDEX(ch_genome_gridss_index_inputs)
            ch_genome_gridss_index = DECOMP_GRIDSS_INDEX.out.extracted_dir

        } else {

            ch_genome_gridss_index = getRefFileChannel('ref_data_genome_gridss_index')

        }
    }

    //
    // Set STAR index , unpack or create if required
    //
    ch_genome_star_index = Channel.empty()
    if (prep_config.require_star_index) {

        if (!params.ref_data_genome_star_index) {

            if(has_alt_contigs)
                error "Refusing to create the STAR index for a ref genome with ALT contigs. Please review https://github.com/alexdobin/STAR docs or contact us on Slack."

            if(!params.ref_data_genome_gtf)
                error "Creating a STAR index requires the appropriate genome transcript annotations as a GTF file. Please contact us on Slack for further information."

            STAR_GENOMEGENERATE(
                ch_genome_fasta,
                file(params.ref_data_genome_gtf),
            )
            ch_genome_star_index = STAR_GENOMEGENERATE.out.index
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        } else if (params.ref_data_genome_star_index.endsWith('.tar.gz')) {

            ch_genome_star_index_inputs = Channel.fromPath(params.ref_data_genome_star_index)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_STAR_INDEX(ch_genome_star_index_inputs)
            ch_genome_star_index = DECOMP_STAR_INDEX.out.extracted_dir

        } else {

            ch_genome_star_index = getRefFileChannel('ref_data_genome_star_index')

        }
    }

    //
    // Set HMF reference data, unpack if required
    //
    ch_hmf_data = Channel.empty()
    if (prep_config.require_hmftools_data) {

        hmf_data_paths = params.hmf_data_paths[params.genome_version]

        if (params.ref_data_hmf_data_path.endsWith('tar.gz')) {

            ch_hmf_data_inputs = Channel.fromPath(params.ref_data_hmf_data_path)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_HMF_DATA(ch_hmf_data_inputs)

            ch_hmf_data = DECOMP_HMF_DATA.out.extracted_dir
                .collect()
                .map { dir_list ->
                    assert dir_list.size() == 1
                    def dirpath = dir_list[0].toUriString()
                    return createDataMap(hmf_data_paths, dirpath)
                }

        } else {

            ch_hmf_data = Channel.value(createDataMap(hmf_data_paths, params.ref_data_hmf_data_path))

        }

        // Set PON paths
        def sequencing_type = RunModes.SequencingType.fromString(params.sequencing_type)

        if(sequencing_type === RunModes.SequencingType.ULTIMA) {

            ch_hmf_data = ch_hmf_data
                .map { d ->
                    if (d.sage_pon_ultima)
                        d.sage_pon = d.sage_pon_ultima

                    if (d.esvee_pon_breakends_ultima)
                        d.esvee_pon_breakends = d.esvee_pon_breakends_ultima

                    if (d.esvee_pon_breakpoints_ultima)
                        d.esvee_pon_breakpoints = d.esvee_pon_breakpoints_ultima

                    return d
                }

        } else if(sequencing_type === RunModes.SequencingType.SBX) {

            ch_hmf_data = ch_hmf_data
                .map { d ->
                    if (d.sage_pon_sbx)
                        d.sage_pon = d.sage_pon_sbx

                    if (d.esvee_pon_breakends_sbx)
                        d.esvee_pon_breakends = d.esvee_pon_breakends_sbx

                    if (d.esvee_pon_breakpoints_sbx)
                        d.esvee_pon_breakpoints = d.esvee_pon_breakpoints_sbx

                    return d
                }

        }

        // Set custom driver gene panel
        if (params.driver_gene_panel) {

            if (pipeline_mode !== RunModes.Pipeline.PANEL_RESOURCE_CREATION) {
                log.info "Using custom driver gene panel: ${params.driver_gene_panel}"
            }

            def custom_driver_panel = file(params.driver_gene_panel, checkIfExists: true)
            ch_hmf_data = ch_hmf_data
                .map { d ->
                    d.driver_gene_panel = custom_driver_panel
                    return d
                }
        }

    }

    //
    // Set panel reference data, unpack if required
    //
    ch_panel_data = Channel.empty()
    if (prep_config.require_panel_data) {

        panel_data_paths_versions = params.panel_data_paths[params.panel]
        panel_data_paths = panel_data_paths_versions[params.genome_version]

        if (params.ref_data_panel_data_path.endsWith('tar.gz')) {

            ch_panel_data_inputs = Channel.fromPath(params.ref_data_panel_data_path)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_PANEL_DATA(ch_panel_data_inputs)

            ch_panel_data = DECOMP_PANEL_DATA.out.extracted_dir
                .collect()
                .map { dir_list ->
                    assert dir_list.size() == 1
                    def dirpath = dir_list[0].toUriString()
                    return createDataMap(panel_data_paths, dirpath)
                }

        } else {

            ch_panel_data = Channel.value(createDataMap(panel_data_paths, params.ref_data_panel_data_path))

        }
    }

    //
    // Write prepared reference data if requested
    //
    if (prepare_reference_only) {

        WRITE_FASTA(ch_genome_fasta)
        WRITE_FAI(ch_genome_fai)
        WRITE_DICT(ch_genome_dict)
        WRITE_IMG(ch_genome_img)
        WRITE_BWA_INDEX(ch_genome_bwamem2_index)
        WRITE_GRIDSS_INDEX(ch_genome_gridss_index)
        WRITE_STAR_INDEX(ch_genome_star_index)

        WRITE_HMF_DATA(ch_hmf_data.map { getDataBaseDirectory(it) })
        WRITE_PANEL_DATA(ch_panel_data.map { getDataBaseDirectory(it) })

    }

    emit:
    genome_fasta         = ch_genome_fasta.first()         // path: genome_fasta
    genome_fai           = ch_genome_fai.first()           // path: genome_fai
    genome_dict          = ch_genome_dict.first()          // path: genome_dict
    genome_img           = ch_genome_img.first()           // path: genome_img
    genome_bwamem2_index = ch_genome_bwamem2_index.first() // path: genome_bwa-mem2_index
    genome_gridss_index  = ch_genome_gridss_index.first()  // path: genome_gridss_index
    genome_star_index    = ch_genome_star_index.first()    // path: genome_star_index
    genome_version       = ch_genome_version               // val:  genome_version

    hmf_data             = ch_hmf_data                     // map:  HMF data paths
    panel_data           = ch_panel_data                   // map:  Panel data paths

    versions             = ch_versions                     // channel: [ versions.yml ]
}

def getConfigForPipelineRun(inputs, pipeline_mode, stages) {

    def has_dna = inputs.any { Inputs.hasTumorDna(it) }
    def has_rna_fastq = inputs.any { Inputs.hasTumorRnaFastq(it) }
    def has_dna_fastq = inputs.any { Inputs.hasTumorDnaFastq(it) || Inputs.hasNormalDnaFastq(it) }

    return [
        require_fasta: true,
        require_fai: true,
        require_dict: true,
        require_img: true,

        require_bwamem2_index: has_dna_fastq && stages.alignment,
        require_star_index: has_rna_fastq && stages.alignment,

        require_gridss_index: has_dna && pipeline_mode === RunModes.Pipeline.WGTS && stages.virusinterpreter,
        require_hmftools_data: true,
        require_panel_data: pipeline_mode === RunModes.Pipeline.TARGETED,
    ]
}

def getConfigForPrepRefOnly(params, log) {

    def ref_data_types = params.ref_data_types
        .tokenize(',')
        .collect { util.Enums.getValidatedEnumFromString(it, RefData.Type) }

    if (
        ref_data_types.contains(RefData.Type.WGS) ||
        ref_data_types.contains(RefData.Type.WTS) ||
        ref_data_types.contains(RefData.Type.TARGETED)
    ) {
        ref_data_types += [
            RefData.Type.FASTA,
            RefData.Type.FAI,
            RefData.Type.DICT,
            RefData.Type.IMG,
            RefData.Type.HMFTOOLS
        ]
    }

    if (ref_data_types.contains(RefData.Type.WGS)) {
        ref_data_types += [RefData.Type.GRIDSS_INDEX]
    }

    if (ref_data_types.contains(RefData.Type.TARGETED)) {
        ref_data_types += [RefData.Type.PANEL]
    }

    def require_fasta = ref_data_types.contains(RefData.Type.FASTA)
    def require_fai = ref_data_types.contains(RefData.Type.FAI)
    def require_dict = ref_data_types.contains(RefData.Type.DICT)
    def require_img = ref_data_types.contains(RefData.Type.IMG)

    def require_bwamem2_index = ref_data_types.contains(RefData.Type.BWAMEM2_INDEX) || ref_data_types.contains(RefData.Type.DNA_ALIGNMENT)
    def require_star_index = ref_data_types.contains(RefData.Type.STAR_INDEX) || ref_data_types.contains(RefData.Type.RNA_ALIGNMENT)

    def require_gridss_index = ref_data_types.contains(RefData.Type.GRIDSS_INDEX)
    def require_hmftools_data = ref_data_types.contains(RefData.Type.HMFTOOLS)
    def require_panel_data = ref_data_types.contains(RefData.Type.PANEL)

    if (require_panel_data) {
        if (params.panel == null) {
            require_panel_data = false
            log.warn "Skipping preparing panel specific reference data as --panel CLI argument was not provided"
        } else if (!RefData.PANELS_DEFINED.contains(params.panel)) {
            require_panel_data = false
            log.warn "Skipping preparing panel specific reference data for custom panel: ${params.panel}"
        }
    }

    return [
        require_fasta: require_fasta,
        require_fai: require_fai,
        require_dict: require_dict,
        require_img: require_img,

        require_bwamem2_index: require_bwamem2_index,
        require_star_index: require_star_index,

        require_gridss_index: require_gridss_index,
        require_hmftools_data: require_hmftools_data,
        require_panel_data: require_panel_data,
    ]
}

def getRefFileChannel(key) {
    def fp = params.get(key) ? file(params.getAt(key)) : []
    return Channel.of(fp)
}

def createDataMap(entries, ref_data_path) {
    return entries
        .collectEntries { name, path ->
            def ref_data_file = path == [] ? [] : getRefdataFile(path, ref_data_path)
            return [name, ref_data_file]
        }
}

def getRefdataFile(filepath, ref_data_path) {
    def data_path_noslash = ref_data_path.toString().replaceAll('/$', '')
    return file("${data_path_noslash}/${filepath}", checkIfExists: true)
}

def getDataBaseDirectory(data) {
    def c = []
    data
        .findAll { it.value }
        .collect { it.value.toUriString().getChars() }
        .transpose()
        .findIndexOf {
            def cs = it.unique()
            if (cs.size() != 1) return true
            c << cs.pop()
            return false
        }
    return file("${c.join('')}")
}
