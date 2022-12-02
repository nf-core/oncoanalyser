//
// This file holds several functions specific to the workflow/oncoanalyser.nf in the nf-core/oncoanalyser pipeline
//

import static groovy.io.FileType.FILES

import nextflow.Channel
import nextflow.Nextflow

import Constants

class WorkflowOncoanalyser {

    //
    // Check and validate parameters
    //
    public static void initialise(params, workflow, log) {

        if (params.genome && params.genomes && ! params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (Constants.genomes_version_37.contains(params.genome)) {
            params.ref_data_genome_version = '37'
        } else if (Constants.genomes_version_38.contains(params.genome)) {
            params.ref_data_genome_version = '38'
        } else if (Constants.genomes_version_38_noalt.contains(params.genome)) {
            params.ref_data_genome_version = '38_noalt'
        } else {
            def genome_version_list_all = Constants.genomes_version_37 + Constants.genomes_version_38
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' is not defined in genome version list.                 \n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${genome_version_list_all.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!params.containsKey('ref_data_hmf_bundle')) {
            if (params.ref_data_genome_version == '37') {
                params.ref_data_hmf_bundle = Constants.hmf_reference_data_37_bundle_path
            } else if (params.ref_data_genome_version == '38') {
                params.ref_data_hmf_bundle = Constants.hmf_reference_data_38_bundle_path
            } else {
                assert false : "Got a bad genome version: ${params.ref_data_genome_version}"
            }
        }

        if (!params.ref_data_genome_fasta) {
            log.error "Genome fasta file not specified with e.g. '--ref_data_genome_fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        // Download region file for collectwgsmetrics if in test mode
        // NOTE(SW): this will be removed as part of the overhaul for testing
        if (workflow.profile.contains('test')) {
            def stage_dir = new File(workflow.workDir.toString(), 'stage/manual/')
            def interval_file = new File(stage_dir, 'collectwgsmetrics.interval_list')
            if (! interval_file.exists()) {
                stage_dir.mkdirs()
                interval_file.createNewFile()
                interval_file << new URL (params.ref_data_wgsmetrics_intervals_url).getText()
            }
            params.ref_data_wgsmetrics_intervals_local = interval_file
        }
    }

    public static prepare_inputs(ch, stub_run, log) {
        return ch
            .splitCsv(header: true, strip: true, sep: '\t')
            .map { [it.id, it] }
            .groupTuple()
            .map { id, inputs ->
                def meta = [
                    'id': id,
                ]
                inputs.each {
                    assert meta.id == it.id
                    // Add subject name if not already present
                    if (meta.containsKey('subject_name')) {
                            assert meta.subject_name == it.subject_name
                    } else {
                            meta.subject_name = it.subject_name
                    }

                    // Set sample name
                    def key = []
                    if (it.sample_type == 'tumor') {
                        key = ['sample_name', 'tumor']
                    } else if (it.sample_type == 'normal') {
                        key = ['sample_name', 'normal']
                    } else {
                        assert false
                    }
                    if (meta.containsKey(key)) {
                        assert meta[key] == it.sample_name
                    } else {
                        meta[key] = it.sample_name
                    }

                    // Add file
                    def key_file = [it.filetype, it.sample_type]
                    assert ! meta.containsKey(key_file)
                    meta[key_file] = it.filepath

                    if (! stub_run) {
                        // For BAM file inputs, require co-located index
                        if (it.filepath.endsWith('.bam')) {
                            def bam_index_fp_str = "${it.filepath}.bai".toString()
                            def bam_index_fp = Nextflow.file(bam_index_fp_str)
                            if (! bam_index_fp.exists()) {
                                log.error "\nERROR: No index found for ${it.filepath}"
                                System.exit(1)
                            }
                        }

                        // For GRIPSS SV VCFs, require co-located index
                        if (it.filetype.startsWith('gripss')) {
                            def vcf_index_fp = Nextflow.file("${it.filepath}.tbi".toString())
                            if (! vcf_index_fp.exists()) {
                                log.error "\nERROR: No index found for ${it.filepath}"
                                System.exit(1)
                            }
                        }
                    }

                    // NOTE(SW): CHECK_SAMPLESHEET curently enforces T/N; this may be relevant in the future
                    //// Set sample type: tumor_normal, tumor_only, normal_only
                    //if (meta.containsKey(['sample_name', 'tumor']) && meta.containsKey(['sample_name', 'normal'])) {
                    //  meta.sample_type = 'tumor_normal'
                    //} else if (meta.containsKey(['sample_name', 'tumor'])) {
                    //  meta.sample_type = 'tumor_only'
                    //} else if (meta.containsKey(['sample_name', 'normal'])) {
                    //  meta.sample_type = 'normal_only'
                    //} else {
                    //  assert false
                    //}
                }

                // For PURPLE only runs, we must get normal sample name from inputs since there is no way to provide this
                // in the samplesheet
                if (meta.containsKey(['cobalt_dir', 'tumor']) && ! meta.containsKey(['sample_name', 'normal'])) {
                    // Discover files
                    def normal_ratio_fps = []
                    new File(meta[['cobalt_dir', 'tumor']])
                        .eachFileMatch(groovy.io.FileType.FILES, ~/.+\.cobalt\.gc\.median\.tsv/, { normal_ratio_fps << it })
                    // Select normal sample file
                    def normal_ratio_fp = normal_ratio_fps
                        .findAll { ! it.getName().contains(meta[['sample_name', 'tumor']]) }
                    assert normal_ratio_fp.size() == 1
                    // Set normal sample name
                    def m = (normal_ratio_fp =~ /.+\/(.+)\.cobalt\.gc\.median\.tsv/)
                    meta[['sample_name', 'normal']] = m[0][1]
                }

                return meta
            }
    }

    public static group_by_meta(Map named_args, ... channels) {
        def r = channels
        // Set position; required to use non-blocking .mix operator
        // NOTE(SW): operating on native list object containing channels
        def i = 0
        r = r
            .collect { ch ->
                def ii = i
                def d = ch.map { data ->
                    def meta = data[0]
                    def values = data[1..-1]
                    return [meta, [position: ii, values: values]]
                }
                i++
                return d
            }

        r = Channel.empty().mix(*r)
        r = r
            .groupTuple(size: channels.size())
            .map { data ->
                def meta = data[0]
                def values_map = data[1]

                def values_list = values_map
                    .sort { it.position }
                    .collect { it.values }
                return [meta, *values_list]
            }

        if (named_args.get('flatten', true)) {
            r = r.map { it.flatten() }
        }
        return r
    }

    // NOTE(SW): function signature required to catch where no named arguments are passed
    public static group_by_meta(... channels) {
        return group_by_meta([:], *channels)
    }

    public static get_input(ch, key) {
        return ch.map { meta -> [meta, meta.getAt(key)] }
    }

    public static restore_meta(ch_output, ch_metas) {
        // NOTE(SW): ch_output must contain a Map in the first position with a key named 'key' that
        // contains the corresponding meta.id value, for example: [val(meta_process), *process_outputs]

        //return Channel.empty()
        //    .mix(
        //        // channel: [val(key), val(meta)]
        //        ch_metas.map { meta -> [meta.id, meta] },
        //        // channel: [val(key), list(process_outputs)]
        //        // NOTE(SW): process meta is discarded
        //        ch_output.map { [it[0].key, it[1..-1]] },
        //    )
        //    .groupTuple(size: 2)
        //    // channel: [meta, *process_outputs]
        //    .map { it.flatten()[1..-1] }

        //def ch_left = ch_metas.map { meta -> [meta.id, meta] }
        //def ch_right = ch_output.map { [it[0].key, it[1..-1]] }
        //return ch_left
        //    .join(ch_right, failOnDuplicate: true)
        //    .map { it[1..-1] }

        def ch_source = ch_metas.map { meta -> [meta.id, meta] }
        def ch_target = ch_output.map { [it[0].key, it[1..-1]] }
        return ch_source
            .cross(ch_target)
            .map { d_meta, d_outputs ->
                def (skey, meta) = d_meta
                def (dkey, outputs) = d_outputs
                return [meta, *outputs]
            }
    }
}
