import SampleSheetFields.FileType
import SampleSheetFields.SampleType
import SampleSheetFields.SequenceType

class Inputs {

    // Files - common
    public static final Map KEY = [

        // Bams
        BAM_DNA_TUMOR:  [FileType.BAM, SampleType.TUMOR, SequenceType.DNA],
        BAM_DNA_NORMAL: [FileType.BAM, SampleType.NORMAL, SequenceType.DNA],
        BAM_DNA_DONOR:  [FileType.BAM, SampleType.DONOR, SequenceType.DNA],
        BAM_RNA_TUMOR:  [FileType.BAM, SampleType.TUMOR, SequenceType.RNA],

        BAI_DNA_TUMOR:  [FileType.BAI, SampleType.TUMOR, SequenceType.DNA],
        BAI_DNA_NORMAL: [FileType.BAI, SampleType.NORMAL, SequenceType.DNA],
        BAI_DNA_DONOR:  [FileType.BAI, SampleType.DONOR, SequenceType.DNA],
        BAI_RNA_TUMOR:  [FileType.BAI, SampleType.TUMOR, SequenceType.RNA],

        // REDUX
        BAM_REDUX_DNA_TUMOR:  [FileType.BAM_REDUX, SampleType.TUMOR, SequenceType.DNA],
        BAM_REDUX_DNA_NORMAL: [FileType.BAM_REDUX, SampleType.NORMAL, SequenceType.DNA],
        BAM_REDUX_DNA_DONOR:  [FileType.BAM_REDUX, SampleType.DONOR, SequenceType.DNA],

        REDUX_TSV_DIR_TUMOR:  [FileType.REDUX_TSV_DIR, SampleType.TUMOR, SequenceType.DNA],
        REDUX_TSV_DIR_NORMAL: [FileType.REDUX_TSV_DIR, SampleType.NORMAL, SequenceType.DNA],
        REDUX_TSV_DIR_DONOR:  [FileType.REDUX_TSV_DIR, SampleType.DONOR, SequenceType.DNA],

        // Other tools
        AMBER_DIR: [FileType.AMBER_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        BAMTOOLS_DIR_TUMOR: [FileType.BAMTOOLS_DIR, SampleType.TUMOR, SequenceType.DNA],
        BAMTOOLS_DIR_NORMAL: [FileType.BAMTOOLS_DIR, SampleType.NORMAL, SequenceType.DNA],

        CHORD_DIR: [FileType.CHORD_DIR, SampleType.TUMOR, SequenceType.DNA],
        COBALT_DIR: [FileType.COBALT_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],
        CUPPA_DIR: [FileType.CUPPA_DIR, SampleType.TUMOR, [SequenceType.DNA, SequenceType.RNA, SequenceType.DNA_RNA]],
        ESVEE_DIR: [FileType.ESVEE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],
        ISOFOX_DIR: [FileType.ISOFOX_DIR, SampleType.TUMOR, SequenceType.RNA],
        LILAC_DIR: [FileType.LILAC_DIR, [SampleType.TUMOR, SampleType.NORMAL, SampleType.TUMOR_NORMAL], [SequenceType.DNA, SequenceType.DNA_RNA]],

        LINX_PLOT_DIR_TUMOR:  [FileType.LINX_PLOT_DIR, SampleType.TUMOR, SequenceType.DNA],
        LINX_ANNO_DIR_TUMOR:  [FileType.LINX_ANNO_DIR, SampleType.TUMOR, SequenceType.DNA],
        LINX_ANNO_DIR_NORMAL: [FileType.LINX_ANNO_DIR, SampleType.NORMAL, SequenceType.DNA],

        SAGE_APPEND_DIR_TUMOR:  [FileType.SAGE_APPEND_DIR, SampleType.TUMOR, SequenceType.DNA_RNA],
        SAGE_APPEND_DIR_NORMAL: [FileType.SAGE_APPEND_DIR, SampleType.NORMAL, SequenceType.DNA_RNA],
        SAGE_DIR_TUMOR:  [FileType.SAGE_DIR, SampleType.TUMOR, SequenceType.DNA],
        SAGE_DIR_NORMAL: [FileType.SAGE_DIR, SampleType.NORMAL, SequenceType.DNA],

        PAVE_DIR_TUMOR:  [FileType.PAVE_DIR, SampleType.TUMOR, SequenceType.DNA],
        PAVE_DIR_NORMAL: [FileType.PAVE_DIR, SampleType.NORMAL, SequenceType.DNA],

        PEACH_DIR: [FileType.PEACH_DIR, SampleType.NORMAL, SequenceType.DNA],
        PURPLE_DIR: [FileType.PURPLE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],
        QSEE_DIR: [FileType.QSEE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],
        SIGS_DIR: [FileType.SIGS_DIR, SampleType.TUMOR, SequenceType.DNA],
        VIRUSINTERPRETER_DIR: [FileType.VIRUSINTERPRETER_DIR, SampleType.TUMOR, SequenceType.DNA],
    ]

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
            case SampleType.TUMOR -> [bam: KEY.BAM_REDUX_DNA_TUMOR, bai: KEY.BAI_DNA_TUMOR]
            case SampleType.NORMAL -> [bam: KEY.BAM_REDUX_DNA_NORMAL, bai: KEY.BAI_DNA_NORMAL]
            case SampleType.DONOR -> [bam: KEY.BAM_REDUX_DNA_DONOR, bai: KEY.BAI_DNA_DONOR]
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
            case SampleType.TUMOR -> KEY.REDUX_TSV_DIR_TUMOR
            case SampleType.NORMAL -> KEY.REDUX_TSV_DIR_NORMAL
            case SampleType.DONOR -> KEY.REDUX_TSV_DIR_DONOR
            default -> throw new IllegalArgumentException("Invalid sample type: ${sample_type}")
        }

        def unwrapped_redux_dir = (redux_dir instanceof List && redux_dir.size() == 1) ? redux_dir[0] : redux_dir

        def selected_redux_dir = preferPipelineOutput(unwrapped_redux_dir, meta, file_key)
        if (!selected_redux_dir)
            return []

        def meta_sample = meta.getOrDefault([sample_type, SequenceType.DNA], [:])
        def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

        // NOTE(LN): toUriString() needs to be called, otherwise nextflow will fail to resolve cloud paths
        def redux_tsvs = nextflow.Nextflow.file("${selected_redux_dir.toUriString()}/${sample_id}.redux.*.tsv*")

        return redux_tsvs
    }


    // Files - SAGE
    private static List resolveSageVcfWithTbi(sage_dir, meta, sample_type) {

        def file_key = switch (sample_type) {
            case SampleType.TUMOR -> KEY.SAGE_DIR_TUMOR
            case SampleType.NORMAL -> KEY.SAGE_DIR_NORMAL
            default -> throw new IllegalArgumentException("Unsupported sample type: ${sample_type}")
        }

        def selected_sage_dir = preferUserProvidedInput(sage_dir, meta, file_key)
        if(!selected_sage_dir)
            return []

        def sample_id = getTumorDnaSampleName(meta)

        def vcf_name = switch (sample_type) {
            case SampleType.TUMOR -> "${sample_id}.sage.somatic.vcf.gz"
            case SampleType.NORMAL -> "${sample_id}.sage.germline.vcf.gz"
            default -> throw new IllegalArgumentException("Unsupported sample type: ${sample_type}")
        }

        return [
            nextflow.Nextflow.file("${selected_sage_dir.toUriString()}/${vcf_name}"),
            nextflow.Nextflow.file("${selected_sage_dir.toUriString()}/${vcf_name}.tbi"),
        ]
    }


    // Files - PURPLE
    public static resolvePurpleSomaticVcf(purple_dir, meta, sample_type) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta, sample_type)}.purple.somatic.vcf.gz")
    }

    public static resolvePurpleSomaticVcf(purple_dir, meta) {
        return resolvePurpleSomaticVcf(purple_dir, meta, 'primary')
    }

    public static resolvePurpleSomaticVcfTbi(purple_dir, meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz.tbi")
    }

    public static resolvePurpleGermlineVcf(purple_dir, meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.germline.vcf.gz")
    }

    public static resolvePurpleSomaticSvVcf(purple_dir, meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.sv.vcf.gz")
    }

    public static resolvePurpleGermlineSvVcf(purple_dir, meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.sv.germline.vcf.gz")
    }


    // Sample records
    public static getTumorDnaSample(meta) { return meta.getOrDefault([SampleType.TUMOR, SequenceType.DNA], [:]) }
    public static getTumorRnaSample(meta) { return meta.getOrDefault([SampleType.TUMOR, SequenceType.RNA], [:]) }
    public static getNormalDnaSample(meta) { return meta.getOrDefault([SampleType.NORMAL, SequenceType.DNA], [:]) }
    public static getDonorDnaSample(meta) { return meta.getOrDefault([SampleType.DONOR, SequenceType.DNA], [:]) }

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

    public static getTumorDnaSampleName(meta) { return getTumorDnaSampleName(meta, 'primary') }

    public static getTumorRnaSampleName(meta) { return getTumorRnaSample(meta)['sample_id'] }
    public static getNormalDnaSampleName(meta) { return getNormalDnaSample(meta)['sample_id'] }
    public static getDonorDnaSampleName(meta) { return getDonorDnaSample(meta)['sample_id'] }


    // Files - Reads/alignments
    public static hasTumorDnaFastq(meta) { return getTumorDnaSample(meta).containsKey(FileType.FASTQ) }
    public static hasTumorDnaBam(meta) { return getTumorDnaSample(meta).containsKey(FileType.BAM) }
    public static hasTumorDnaReduxBam(meta) { return getTumorDnaSample(meta).containsKey(FileType.BAM_REDUX) }

    public static hasNormalDnaFastq(meta) { return getNormalDnaSample(meta).containsKey(FileType.FASTQ) }
    public static hasNormalDnaBam(meta) { return getNormalDnaSample(meta).containsKey(FileType.BAM) }
    public static hasNormalDnaReduxBam(meta) { return getNormalDnaSample(meta).containsKey(FileType.BAM_REDUX) }

    public static hasDonorDnaFastq(meta) { return getDonorDnaSample(meta).containsKey(FileType.FASTQ) }
    public static hasDonorDnaBam(meta) { return getDonorDnaSample(meta).containsKey(FileType.BAM) }
    public static hasDonorDnaReduxBam(meta) { return getDonorDnaSample(meta).containsKey(FileType.BAM_REDUX) }

    public static hasTumorRnaFastq(meta) { return getTumorRnaSample(meta).containsKey(FileType.FASTQ) }
    public static hasTumorRnaBam(meta) { return getTumorRnaSample(meta).containsKey(FileType.BAM) }


    // Status
    public static hasTumorDna(meta) { return hasTumorDnaBam(meta) || hasTumorDnaReduxBam(meta) || hasTumorDnaFastq(meta) }
    public static hasNormalDna(meta) { return hasNormalDnaBam(meta) || hasNormalDnaReduxBam(meta) || hasNormalDnaFastq(meta) }
    public static hasDonorDna(meta) { return hasDonorDnaBam(meta) || hasDonorDnaReduxBam(meta) || hasDonorDnaFastq(meta) }
    public static hasTumorRna(meta) { return hasTumorRnaBam(meta) || hasTumorRnaFastq(meta) }
}
