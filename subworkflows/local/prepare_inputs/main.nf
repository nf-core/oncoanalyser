//
// Prepare inputs (tests only)
//

// NOTE(SW): inputs for the pipeline are prepared outside of NF
// workflow/channels to allow higher-level conditionals, however nf-test
// well-formed meta (including Constants) that can only be made available
// through running workflows/processes with 'setup'. Hence, this subworkflow
// isn't used in the main pipeline and is only used for execution of tests.

import Utils

workflow PREPARE_INPUTS {
    take:
    input_fp_str

    main:
    ch_inputs = Channel.fromList(
        Utils.parseInput(input_fp_str, workflow.stubRun, log)
    )

    emit:
    inputs = ch_inputs // channel: [ meta ]
}
