// NOTE: Full file not provided here because only one targeted block needs change and
// minimizing file churn is best. Replace ONLY the block below in your existing file.

// Select input sources and sort for fusion annotation
// channel: runnable: [ meta, neo_finder_dir, tumor_bam, tumor_bai ]
// channel: skip: [ meta ]
ch_isofox_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
    ch_finder_out,
    ch_tumor_rna_bam,
)
    .map { meta, neo_finder_dir, tumor_bam, tumor_bai ->
        return [
            meta,
            neo_finder_dir,
            params.realign_bam
                ? tumor_bam
                : Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
            params.realign_bam
                ? tumor_bai
                : Utils.selectCurrentOrExisting(tumor_bai, meta, Constants.INPUT.BAI_RNA_TUMOR),
        ]
    }
    .branch { meta, neo_finder_dir, tumor_bam, tumor_bai ->
        runnable: Utils.hasTumorRna(meta)
            return [meta, neo_finder_dir, tumor_bam, tumor_bai]
        skip: true
            return meta
    }
