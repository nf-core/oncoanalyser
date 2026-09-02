//
// WISP estimates tumor purity in longitudinal samples using WGS data of the primary
//

include { WISP } from '../../../modules/local/wisp/main'

workflow WISP_ANALYSIS {
    take:
    // Sample data
    ch_inputs                  // channel: [mandatory] [ meta ]
    ch_redux_dir               // channel: [mandatory] [ meta, redux_dir ]
    ch_amber_dir               // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_dir              // channel: [mandatory] [ meta, cobalt_dir ]
    ch_sage_append_dir_somatic // channel: [mandatory] [ meta, sage_append_dir ]

    // Reference data
    genome_fasta               // channel: [mandatory] /path/to/genome_fasta
    genome_fai                 // channel: [mandatory] /path/to/genome_fai

    // Params
    targeted_mode              // boolean: [mandatory] Set targeted mode

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir (primary), amber_dir (primary), normal_aln (primary), redux_dir (longitudinal), amber_dir (longitudinal), cobalt_dir (longitudinal), sage_append_dir (longitudinal) ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir,
        ch_amber_dir,
        ch_cobalt_dir,
        ch_sage_append_dir_somatic,
    )
        .map { meta, longitudinal_redux_dir, longitudinal_amber_dir, longitudinal_cobalt_dir, longitudinal_sage_append_dir ->

            def primary_normal_redux_dir = Utils.getInput(meta, Constants.INPUT.REDUX_DIR_NORMAL)
            def (primary_normal_aln, _primary_normal_idx) = Utils.getNormalReduxDirAlignment(meta, primary_normal_redux_dir)

            def primary_purple_dir = Utils.getInput(meta, Constants.INPUT.PURPLE_DIR)
            def primary_amber_dir = Utils.getInput(meta, Constants.INPUT.AMBER_DIR)

            def longitudinal_redux_dir_selected = Utils.selectCurrentOrExisting(longitudinal_redux_dir, meta, Constants.INPUT.REDUX_DIR_TUMOR)

            return [
              meta,
              primary_purple_dir,
              primary_amber_dir,
              primary_normal_aln,
              longitudinal_redux_dir_selected,
              longitudinal_amber_dir,
              longitudinal_cobalt_dir,
              longitudinal_sage_append_dir,
            ]
        }
        .branch { meta, primary_purple_dir, primary_amber_dir, primary_normal_aln, longitudinal_redux_dir, longitudinal_amber_dir, longitudinal_cobalt_dir, longitudinal_sage_append_dir ->

            def runnable
            if (targeted_mode) {
                runnable = primary_purple_dir && longitudinal_sage_append_dir
            } else {
                runnable = primary_purple_dir && longitudinal_sage_append_dir && longitudinal_cobalt_dir
            }

            runnable: runnable
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_wisp, purple_dir (primary), amber_dir (primary), normal_aln (primary), redux_dir (longitudinal), amber_dir (longitudinal), cobalt_dir (longitudinal), sage_append_dir (longitudinal) ]
    ch_wisp_inputs = ch_inputs_sorted.runnable
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def meta_wisp = [
                key: meta.group_id,
                id: meta.group_id,
                subject_id: meta.subject_id,
                primary_id: Utils.getTumorDnaSampleName(meta, primary: true),
                longitudinal_id: Utils.getTumorDnaSampleName(meta, primary: false),
            ]

            return [meta_wisp] + inputs
        }

    // Run process
    WISP(
        ch_wisp_inputs,
        genome_fasta,
        genome_fai,
        targeted_mode,
    )
}
