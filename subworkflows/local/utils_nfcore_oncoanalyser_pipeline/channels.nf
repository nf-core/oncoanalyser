//
// Input channel builders for the nf-core/oncoanalyser pipeline
//

include { getDonorDnaSamples  } from './accessors'
include { getNormalDnaSample  } from './accessors'
include { getTumorDnaSample   } from './accessors'
include { getTumorRnaSample   } from './accessors'
include { hasDonorDnaAlns     } from './accessors'
include { hasDonorDnaFastqs   } from './accessors'
include { hasInput            } from './accessors'
include { hasNormalDnaAln     } from './accessors'
include { hasNormalDnaFastq   } from './accessors'
include { hasTumorDnaAln      } from './accessors'
include { hasTumorDnaFastq    } from './accessors'
include { hasTumorRnaFastq    } from './accessors'
include { FileType            } from './types'

def getDnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ case_record ]
    def ch_inputs_tumor_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasTumorDnaAln(case_record)
            runnable: hasTumorDnaFastq(case_record) && ! has_existing
            skip: true
        }

    def ch_inputs_normal_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasNormalDnaAln(case_record)
            runnable: hasNormalDnaFastq(case_record) && ! has_existing
            skip: true
        }

    def ch_inputs_donor_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasDonorDnaAlns(case_record)
            runnable: hasDonorDnaFastqs(case_record) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ case_record, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { case_record -> [case_record, getTumorDnaSample(case_record), 'tumor'] },
            ch_inputs_normal_sorted.runnable.map { case_record -> [case_record, getNormalDnaSample(case_record), 'normal'] },
            ch_inputs_donor_sorted.runnable.flatMap { case_record -> getDonorDnaSamples(case_record).collect { donor -> [case_record, donor, 'donor'] } },
        )
        .flatMap { case_record, sample, sample_type ->
            sample.files.getAt(FileType.FASTQ)
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

                    return [case_record, fastq_info, fastq.read_fwd, fastq.read_rev ?: null]
                }
        }

    return channel.empty()
        .mix(
            ch_fastqs,
            ch_inputs_tumor_sorted.skip.map { case_record -> [case_record, [:], null, null] },
            ch_inputs_normal_sorted.skip.map { case_record -> [case_record, [:], null, null] },
            ch_inputs_donor_sorted.skip.map { case_record -> [case_record, [:], null, null] },
        )
}

def getRnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ case_record ]
    def ch_inputs_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasInput(getTumorRnaSample(case_record), FileType.ALN)
            runnable: hasTumorRnaFastq(case_record) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ case_record, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = ch_inputs_sorted.runnable
        .flatMap { case_record ->
            def sample = getTumorRnaSample(case_record)
            sample.files
                .getAt(FileType.FASTQ)
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
            ch_inputs_sorted.skip.map { case_record -> [case_record, [:], null, null] },
        )
}
