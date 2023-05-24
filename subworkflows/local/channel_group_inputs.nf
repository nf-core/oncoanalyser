import Constants
import Utils

workflow CHANNEL_GROUP_INPUTS {
    take:
        ch_inputs

    main:
        // channel (present): [val(meta)]
        // channel (absent): [val(meta)]
        ch_inputs_wgs = ch_inputs
            .branch { meta ->
                def key_tumor = [Constants.FileType.BAM, Constants.SampleType.TUMOR, Constants.SequenceType.WGS]
                def key_normal = [Constants.FileType.BAM, Constants.SampleType.NORMAL, Constants.SequenceType.WGS]
                present: meta.containsKey(key_tumor) && meta.containsKey(key_normal)
                    return meta
                absent: true
                    return meta
            }

        // channel (present): [val(meta)]
        // channel (absent): [val(meta)]
        ch_inputs_wts = ch_inputs
            .branch { meta ->
                def key = [Constants.FileType.BAM, Constants.SampleType.TUMOR, Constants.SequenceType.WTS]
                present: meta.containsKey(key)
                    return meta
                absent: ! meta.containsKey(key)
                    return meta
            }

    emit:
        wgs_present = ch_inputs_wgs.present
        wgs_absent = ch_inputs_wgs.absent

        wts_present = ch_inputs_wts.present
        wts_absent = ch_inputs_wts.absent
}
