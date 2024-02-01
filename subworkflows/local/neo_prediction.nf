//
// XXX
//

import Constants
import Utils

include { ISOFOX               } from '../../modules/local/isofox/main'
include { LILAC                } from '../../modules/local/lilac/main'
include { NEO as NEO_PREDICTOR } from '../../modules/local/neo/main'
include { NEO as NEO_SCORER    } from '../../modules/local/neo/main'

workflow NEO_PREDICTION{
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_isofox              // channel: [mandatory] [ meta, isofox_dir ]
        ch_purple              // channel: [mandatory] [ meta, purple_dir ]
        ch_sage_somatic_append // channel: [mandatory] [ meta, sage_append_vcf ]
        ch_lilac               // channel: [mandatory] [ meta, lilac_dir ]
        ch_linx                // channel: [mandatory] [ meta, linx_dir ]

        //// Reference data
        //genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        //genome_fai             // channel: [mandatory] /path/to/genome_fai
        //genome_dict            // channel: [mandatory] /path/to/genome_dict

        // other reference data, or placeholders

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Neo prediction
        // 1. select input sources after combine required channels (i.e ch_purple, ch_linx)
        // 2. get runnable subjects/inputs
        // 3. format input channel
        // 4. run process NEO_PREDICTOR
        // 5. restore meta, set skip entries

        // Feeding the Neo process raw inputs for demo purposes only
        NEO_PREDICTOR(ch_inputs)

        // Isofox annotation
        // 1. take outputs from NEO_PREDICTOR
        // 2. format input channel
        // 3. run process ISOFOX (also requires changes to process to enable alt. run mode)
        // 4. restore meta

        // Assuming SAGE append will be appropriately handled upstream. Need to discuss some options:
        //   * -bqr_enabled true
        //   * -max_read_depth 100000

        // Neo score
        // 1. select input sources after combine required channels (i.e ch_purple, ch_lilac, ch_isofox, above process outputs)
        // 2. get runnable subjects/inputs
        // 3. format input channel
        // 4. run NEO_SCORER

        // Feeding the Neo process raw inputs for demo purposes only
        NEO_SCORER(ch_inputs)

    emit:
        versions = ch_versions // channel: [ versions.yml ]
}
