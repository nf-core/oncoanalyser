class Inputs {

    // Sample records
    public static getTumorDnaSample(meta) {
        return meta.getOrDefault([SampleMeta.SampleType.TUMOR, SampleMeta.SequenceType.DNA], [:])
    }

    public static getTumorRnaSample(meta) {
        return meta.getOrDefault([SampleMeta.SampleType.TUMOR, SampleMeta.SequenceType.RNA], [:])
    }

    public static getNormalDnaSample(meta) {
        return meta.getOrDefault([SampleMeta.SampleType.NORMAL, SampleMeta.SequenceType.DNA], [:])
    }

    public static getDonorDnaSample(meta) {
        return meta.getOrDefault([SampleMeta.SampleType.DONOR, SampleMeta.SequenceType.DNA], [:])
    }

    // Sample names
    public static getTumorDnaSampleName(meta, sample_type) {
        def meta_sample = getTumorDnaSample(meta)
        def sample_id

        // NOTE(LN): Sample type is a string (and not a boolean) so that it is obvious what sample type is being
        // retrieved when this method is called. This wouldn't be necessary if Nextflow had IDE support, where
        // the argument name would be displayed as an overlay
        if(sample_type == 'primary') {
            sample_id = meta_sample['sample_id']
        } else if (sample_type == 'longitudinal') {
            sample_id = meta_sample['longitudinal_sample_id']
        } else {
            throw new IllegalArgumentException("`sample_type` must be 'primary' or 'longitudinal'")
        }

        return sample_id
    }

    public static getTumorDnaSampleName(meta) {
        return getTumorDnaSampleName(meta, 'primary')
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
        return getTumorDnaSample(meta).getOrDefault(SampleMeta.FileType.FASTQ, null)
    }

    public static getTumorDnaBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(SampleMeta.FileType.BAM, null)
    }

    public static getTumorDnaReduxBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(SampleMeta.FileType.BAM_REDUX, null)
    }


    public static hasTumorDnaFastq(meta) {
        return getTumorDnaFastq(meta) !== null
    }

    public static hasTumorDnaBam(meta) {
        return getTumorDnaBam(meta) !== null
    }

    public static hasTumorDnaReduxBam(meta) {
        return getTumorDnaReduxBam(meta) !== null
    }


    // Files - Normal DNA
    public static getNormalDnaFastq(meta) {
        return getNormalDnaSample(meta).getOrDefault(SampleMeta.FileType.FASTQ, null)
    }

    public static getNormalDnaBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(SampleMeta.FileType.BAM, null)
    }

    public static getNormalDnaReduxBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(SampleMeta.FileType.BAM_REDUX, null)
    }


    public static hasNormalDnaFastq(meta) {
        return getNormalDnaFastq(meta) !== null
    }

    public static hasNormalDnaBam(meta) {
        return getNormalDnaBam(meta) !== null
    }

    public static hasNormalDnaReduxBam(meta) {
        return getNormalDnaReduxBam(meta) !== null
    }

    // Files - Donor DNA
    public static getDonorDnaFastq(meta) {
        return getDonorDnaSample(meta).getOrDefault(SampleMeta.FileType.FASTQ, null)
    }

    public static getDonorDnaBam(meta) {
        return getDonorDnaSample(meta).getOrDefault(SampleMeta.FileType.BAM, null)
    }

    public static getDonorDnaReduxBam(meta) {
        return getDonorDnaSample(meta).getOrDefault(SampleMeta.FileType.BAM_REDUX, null)
    }


    public static hasDonorDnaFastq(meta) {
        return getDonorDnaFastq(meta) !== null
    }

    public static hasDonorDnaBam(meta) {
        return getDonorDnaBam(meta) !== null
    }

    public static hasDonorDnaReduxBam(meta) {
        return getDonorDnaReduxBam(meta) !== null
    }


    // Files - Tumor RNA
    public static getTumorRnaFastq(meta) {
        return getTumorRnaSample(meta).getOrDefault(SampleMeta.FileType.FASTQ, null)
    }

    public static getTumorRnaBam(meta) {
        return getTumorRnaSample(meta).getOrDefault(SampleMeta.FileType.BAM, null)
    }

    public static hasTumorRnaFastq(meta) {
        return getTumorRnaFastq(meta) !== null
    }

    public static hasTumorRnaBam(meta) {
        return getTumorRnaBam(meta) !== null
    }


    // Status
    public static hasTumorDna(meta) {
        return hasTumorDnaBam(meta) || hasTumorDnaReduxBam(meta) || hasTumorDnaFastq(meta)
    }

    public static hasNormalDna(meta) {
        return hasNormalDnaBam(meta) || hasNormalDnaReduxBam(meta) || hasNormalDnaFastq(meta)
    }

    public static hasDonorDna(meta) {
        return hasDonorDnaBam(meta) || hasDonorDnaReduxBam(meta) || hasDonorDnaFastq(meta)
    }

    public static hasTumorRna(meta) {
        return hasTumorRnaBam(meta) || hasTumorRnaFastq(meta)
    }

    // Files - PURPLE
    //
    // NOTE(LN): We construct the output file paths based on the tool output dir. This so that the user can resume a run by providing
    // purple_dir in the sample sheet.
    //
    // Processes in theory could emit these output files, and downstream processes could consume the output file channels. However,
    // this would make resuming messy as multiple e.g. PURPLE output files would need to be provided in the sample sheet

    public static getPurpleSomaticVcf(meta, purple_dir, sample_type) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta, sample_type)}.purple.somatic.vcf.gz")
    }

    public static getPurpleSomaticVcf(meta, purple_dir) {
        return getPurpleSomaticVcf(meta, purple_dir, 'primary')
    }

    public static getPurpleSomaticVcfTbi(meta, purple_dir) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz.tbi")
    }

    public static getPurpleGermlineVcf(meta, purple_dir) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.germline.vcf.gz")
    }

    public static getPurpleSomaticSvVcf(meta, purple_dir) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.sv.vcf.gz")
    }

    public static getPurpleGermlineSvVcf(meta, purple_dir) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.sv.germline.vcf.gz")
    }

    // Files - REDUX
    public static resolveReduxBamBai(redux_bam_bai, meta, sample_type) {

        def key_map = [
            (SampleMeta.SampleType.TUMOR): [
                bam: SampleMeta.INPUT.BAM_REDUX_DNA_TUMOR,
                bai: SampleMeta.INPUT.BAI_DNA_TUMOR,
            ],

            (SampleMeta.SampleType.NORMAL): [
                bam: SampleMeta.INPUT.BAM_REDUX_DNA_NORMAL,
                bai: SampleMeta.INPUT.BAI_DNA_NORMAL,
            ],

            (SampleMeta.SampleType.DONOR): [
                bam: SampleMeta.INPUT.BAM_REDUX_DNA_DONOR,
                bai: SampleMeta.INPUT.BAI_DNA_DONOR,
            ],
        ]

        def keys = key_map[sample_type]

        def (bam, bai) = redux_bam_bai

        return [
            preferUserProvidedInput(bam, meta, keys.bam),
            preferPipelineOutput(bai, meta, keys.bai),
        ]
    }

    public static resolveReduxTsvFiles(redux_tsvs, meta, sample_type) {

        def key_map = [
            (SampleMeta.SampleType.TUMOR): [
                bqr_tsv: SampleMeta.INPUT.REDUX_BQR_TSV_TUMOR,
                dup_freq_tsv: SampleMeta.INPUT.REDUX_DUP_FREQ_TSV_TUMOR,
                jitter_tsv: SampleMeta.INPUT.REDUX_JITTER_TSV_TUMOR,
                ms_tsv: SampleMeta.INPUT.REDUX_MS_TSV_TUMOR,
            ],

            (SampleMeta.SampleType.NORMAL): [
                bqr_tsv: SampleMeta.INPUT.REDUX_BQR_TSV_NORMAL,
                dup_freq_tsv: SampleMeta.INPUT.REDUX_DUP_FREQ_TSV_NORMAL,
                jitter_tsv: SampleMeta.INPUT.REDUX_JITTER_TSV_NORMAL,
                ms_tsv: SampleMeta.INPUT.REDUX_MS_TSV_NORMAL,
            ],

            (SampleMeta.SampleType.DONOR): [
                bqr_tsv: SampleMeta.INPUT.REDUX_BQR_TSV_DONOR,
                dup_freq_tsv: SampleMeta.INPUT.REDUX_DUP_FREQ_TSV_DONOR,
                jitter_tsv: SampleMeta.INPUT.REDUX_JITTER_TSV_DONOR,
                ms_tsv: SampleMeta.INPUT.REDUX_MS_TSV_DONOR,
            ],
        ]

        def keys = key_map[sample_type]

        def (bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv) = redux_tsvs

        return [
            preferPipelineOutput(bqr_tsv, meta, keys.bqr_tsv),
            preferPipelineOutput(dup_freq_tsv, meta, keys.dup_freq_tsv),
            preferPipelineOutput(jitter_tsv, meta, keys.jitter_tsv),
            preferPipelineOutput(ms_tsv, meta, keys.ms_tsv),
        ].findAll { it != [] }
    }

    // Misc
    public static getInput(meta, key) {

        def result = []
        def (file_type, sample_types, sequence_type) = key

        for (key_sample in [sample_types, sequence_type].combinations()) {
            if (meta.containsKey(key_sample) && meta[key_sample].containsKey(file_type)) {
                result = meta[key_sample].get(file_type)
                break
            }
        }
        return result
    }

    public static hasExistingInput(meta, key) {
        return getInput(meta, key) != []
    }

    public static preferUserProvidedInput(pipeline_path, meta, key) {
        // Allows the pipeline to start from downstream steps, e.g. running ORANGE from existing pipeline outputs
        def user_provided_path = getInput(meta, key)
        return user_provided_path ?: pipeline_path
    }

    public static preferPipelineOutput(pipeline_path, meta, key) {
        def user_provided_path = getInput(meta, key)
        return pipeline_path ?: user_provided_path
    }
}
