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

    // Files - common
    public static getInput(meta, key) {

        def result = []
        def (file_type, sample_types, sequence_type) = key

        for (sample_key in [sample_types, sequence_type].combinations()) {
            if (meta.containsKey(sample_key) && meta[sample_key].containsKey(file_type)) {
                result = meta[sample_key].get(file_type)
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

    // Files - REDUX
    public static List resolveReduxBamBai(redux_bam_bai, meta, sample_type) {

        def file_keys = switch (sample_type) {
            case SampleMeta.SampleType.TUMOR -> [bam: SampleMeta.INPUT.BAM_REDUX_DNA_TUMOR, bai: SampleMeta.INPUT.BAI_DNA_TUMOR]
            case SampleMeta.SampleType.NORMAL -> [bam: SampleMeta.INPUT.BAM_REDUX_DNA_NORMAL, bai: SampleMeta.INPUT.BAI_DNA_NORMAL]
            case SampleMeta.SampleType.DONOR -> [bam: SampleMeta.INPUT.BAM_REDUX_DNA_DONOR, bai: SampleMeta.INPUT.BAI_DNA_DONOR]
            default -> throw new IllegalArgumentException("Invalid sample type: ${sample_type}")
        }

        def (bam, bai) = redux_bam_bai

        return [
            preferUserProvidedInput(bam, meta, file_keys.bam),
            preferPipelineOutput(bai, meta, file_keys.bai),
        ]
    }

    public static List resolveReduxTsvFiles(redux_dir, meta, sample_type) {

        // NOTE(LN): Get the REDUX TSV files as a glob
        //
        // This avoids passing REDUX dir (containing the BAM and TSV files) to downstream processes because in
        // cloud environments, this would mean the BAMs are copied to the VM running that downstream process.
        // When the BAMs are not needed as input, this would result in the VM requiring more disk space than necessary.

        def file_key = switch (sample_type) {
            case SampleMeta.SampleType.TUMOR -> SampleMeta.INPUT.REDUX_TSV_DIR_TUMOR
            case SampleMeta.SampleType.NORMAL -> SampleMeta.INPUT.REDUX_TSV_DIR_NORMAL
            case SampleMeta.SampleType.DONOR -> SampleMeta.INPUT.REDUX_TSV_DIR_DONOR
            default -> throw new IllegalArgumentException("Invalid sample type: ${sample_type}")
        }

        def unwrapped_redux_dir = (redux_dir instanceof List && redux_dir.size() == 1) ? redux_dir[0] : redux_dir

        def selected_redux_dir = preferPipelineOutput(unwrapped_redux_dir, meta, file_key)
        if (!selected_redux_dir)
            return []

        def meta_sample = meta.getOrDefault([sample_type, SampleMeta.SequenceType.DNA], [:])
        def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

        def redux_tsvs = nextflow.Nextflow.files("${selected_redux_dir}/${sample_id}.redux.*.tsv*")
        if (!redux_tsvs)
            return []

        return redux_tsvs
    }

    // Files - SAGE
    private static List resolveSageVcfWithTbi(sage_dir, meta, sample_type) {

        def file_key = switch (sample_type) {
            case SampleMeta.SampleType.TUMOR -> SampleMeta.INPUT.SAGE_DIR_TUMOR
            case SampleMeta.SampleType.NORMAL -> SampleMeta.INPUT.SAGE_DIR_NORMAL
            default -> throw new IllegalArgumentException("Unsupported sample type: ${sample_type}")
        }

        def selected_sage_dir = preferUserProvidedInput(sage_dir, meta, file_key)
        if(!selected_sage_dir)
            return []

        def sage_dir_path = nextflow.Nextflow.file(selected_sage_dir)

        def sample_id = getTumorDnaSampleName(meta)

        def vcf_name = switch (sample_type) {
            case SampleMeta.SampleType.TUMOR -> "${sample_id}.sage.somatic.vcf.gz"
            case SampleMeta.SampleType.NORMAL -> "${sample_id}.sage.germline.vcf.gz"
            default -> throw new IllegalArgumentException("Unsupported sample type: ${sample_type}")
        }

        return [ sage_dir_path.resolve(vcf_name), sage_dir_path.resolve("${vcf_name}.tbi") ]
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
}
