package sample

import samplesheet.FileType
import samplesheet.SampleType
import samplesheet.SequenceType

import java.nio.file.Path

class Inputs {

    public static Object get(Map meta, FileKey key) {

        def result = []

        def file_type = key.fileType
        def sample_types = key.sampleTypes
        def sequence_types = key.sequenceTypes

        def sample_sequence_type_combinations = [sample_types, sequence_types].combinations()

        for (sample_key in sample_sequence_type_combinations) {
            if (meta.containsKey(sample_key) && meta[sample_key].containsKey(file_type)) {
                result = meta[sample_key].get(file_type)
                break
            }
        }
        return result
    }

    public static boolean hasExisting(Map meta, FileKey key) {
        return get(meta, key) != []
    }

    public static Object preferUserProvidedInput(Object pipeline_path, Map meta, FileKey key) {
        // Allows the pipeline to start from downstream steps, e.g. running ORANGE from existing pipeline outputs
        def user_provided_path = get(meta, key)
        return user_provided_path ?: pipeline_path
    }

    public static Object preferPipelineOutput(Object pipeline_path, Map meta, FileKey key) {
        def user_provided_path = get(meta, key)
        return pipeline_path ?: user_provided_path
    }


    // Files - REDUX
    public static List<Object> resolveReduxBamBai(List<Object> redux_bam_bai, Map meta, SampleType sample_type) {

        def file_keys = switch (sample_type) {
            case SampleType.TUMOR -> [bam: FileKey.BAM_REDUX_DNA_TUMOR, bai: FileKey.BAI_DNA_TUMOR]
            case SampleType.NORMAL -> [bam: FileKey.BAM_REDUX_DNA_NORMAL, bai: FileKey.BAI_DNA_NORMAL]
            case SampleType.DONOR -> [bam: FileKey.BAM_REDUX_DNA_DONOR, bai: FileKey.BAI_DNA_DONOR]
            default -> throw new IllegalArgumentException("Invalid sample type: ${sample_type}")
        }

        def (bam, bai) = redux_bam_bai

        return [
            preferUserProvidedInput(bam, meta, file_keys.bam),
            preferPipelineOutput(bai, meta, file_keys.bai),
        ]
    }

    public static List<Path> resolveReduxTsvFiles(Object redux_dir, Map meta, SampleType sample_type) {

        // NOTE(LN): This avoids passing REDUX dir (containing the BAM and TSV files) to downstream processes because in
        // cloud environments, this would mean the BAMs are copied to the VM running that downstream process.
        // When the BAMs are not needed as input, this would result in the VM requiring more disk space than necessary.

        def file_key = switch (sample_type) {
            case SampleType.TUMOR -> FileKey.REDUX_TSV_DIR_TUMOR
            case SampleType.NORMAL -> FileKey.REDUX_TSV_DIR_NORMAL
            case SampleType.DONOR -> FileKey.REDUX_TSV_DIR_DONOR
            default -> throw new IllegalArgumentException("Invalid sample type: ${sample_type}")
        }

        def unwrapped_redux_dir = (redux_dir instanceof List && redux_dir.size() == 1) ? redux_dir[0] : redux_dir

        def selected_redux_dir = preferPipelineOutput(unwrapped_redux_dir, meta, file_key)
        if (!selected_redux_dir)
            return []

        def meta_sample = meta.getOrDefault([sample_type, SequenceType.DNA], [:])
        def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

        // NOTE(LN):
        // -  Ideally we would use nextflow.Nextflow.file(...).resolve(...) but that returns a Path with wildcard paths
        //   NOT expanded. However, we want List<Path> where wildcard paths are expanded.
        // - toUriString() needs to be called, otherwise nextflow will fail to resolve cloud paths.
        def redux_dir_uri = nextflow.Nextflow.file(selected_redux_dir).toUriString()
        def redux_tsvs = nextflow.Nextflow.file("${redux_dir_uri}/${sample_id}.redux.*.tsv*")

        return redux_tsvs
    }


    // Files - SAGE
    private static List<Path> resolveSageVcfWithTbi(Object sage_dir, Map meta, SampleType sample_type) {

        def file_key = switch (sample_type) {
            case SampleType.TUMOR -> FileKey.SAGE_DIR_TUMOR
            case SampleType.NORMAL -> FileKey.SAGE_DIR_NORMAL
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
            nextflow.Nextflow.file(selected_sage_dir).resolve(vcf_name),
            nextflow.Nextflow.file(selected_sage_dir).resolve(vcf_name + '.tbi'),
        ]
    }


    // Files - PURPLE
    public static Path resolvePurpleSomaticVcf(Path purple_dir, Map meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz")
    }

    public static Path resolvePurpleSomaticVcfTbi(Path purple_dir, Map meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz.tbi")
    }

    public static Path resolvePurpleGermlineVcf(Path purple_dir, Map meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.germline.vcf.gz")
    }

    public static Path resolvePurpleSomaticSvVcf(Path purple_dir, Map meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.sv.vcf.gz")
    }

    public static Path resolvePurpleGermlineSvVcf(Path purple_dir, Map meta) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.sv.germline.vcf.gz")
    }


    // Sample records
    public static Map getTumorDnaSample(Map meta) { return meta.getOrDefault([SampleType.TUMOR, SequenceType.DNA], [:]) }
    public static Map getTumorRnaSample(Map meta) { return meta.getOrDefault([SampleType.TUMOR, SequenceType.RNA], [:]) }
    public static Map getNormalDnaSample(Map meta) { return meta.getOrDefault([SampleType.NORMAL, SequenceType.DNA], [:]) }
    public static Map getDonorDnaSample(Map meta) { return meta.getOrDefault([SampleType.DONOR, SequenceType.DNA], [:]) }

    // Sample names
    public static String getTumorDnaSampleName(Map meta) { return getTumorDnaSample(meta)['sample_id'] }
    public static String getTumorDnaSampleNamePrimary(Map meta) { return getTumorDnaSampleName(meta) } // NOTE(LN): Alias method used for clarity
    public static String getTumorDnaSampleNameLongitudinal(Map meta) { return getTumorDnaSample(meta)['longitudinal_sample_id'] }

    public static String getTumorRnaSampleName(Map meta) { return getTumorRnaSample(meta)['sample_id'] }
    public static String getNormalDnaSampleName(Map meta) { return getNormalDnaSample(meta)['sample_id'] }
    public static String getDonorDnaSampleName(Map meta) { return getDonorDnaSample(meta)['sample_id'] }

    public static String getTumorRnaSampleOutputId(Map meta) {
        // For isofox output files, use the tumor DNA sample ID, but fall back to the tumor RNA sample ID when running in RNA-only mode
        def tumor_dna_id = getTumorDnaSampleName(meta)
        def tumor_rna_id = getTumorRnaSampleName(meta)
        return (tumor_dna_id) ? tumor_dna_id : tumor_rna_id
    }


    // Files - Reads/alignments
    public static boolean hasTumorDnaFastq(Map meta) { return getTumorDnaSample(meta).containsKey(FileType.FASTQ) }
    public static boolean hasTumorDnaBam(Map meta) { return getTumorDnaSample(meta).containsKey(FileType.BAM) }
    public static boolean hasTumorDnaReduxBam(Map meta) { return getTumorDnaSample(meta).containsKey(FileType.BAM_REDUX) }

    public static boolean hasNormalDnaFastq(Map meta) { return getNormalDnaSample(meta).containsKey(FileType.FASTQ) }
    public static boolean hasNormalDnaBam(Map meta) { return getNormalDnaSample(meta).containsKey(FileType.BAM) }
    public static boolean hasNormalDnaReduxBam(Map meta) { return getNormalDnaSample(meta).containsKey(FileType.BAM_REDUX) }

    public static boolean hasDonorDnaFastq(Map meta) { return getDonorDnaSample(meta).containsKey(FileType.FASTQ) }
    public static boolean hasDonorDnaBam(Map meta) { return getDonorDnaSample(meta).containsKey(FileType.BAM) }
    public static boolean hasDonorDnaReduxBam(Map meta) { return getDonorDnaSample(meta).containsKey(FileType.BAM_REDUX) }

    public static boolean hasTumorRnaFastq(Map meta) { return getTumorRnaSample(meta).containsKey(FileType.FASTQ) }
    public static boolean hasTumorRnaBam(Map meta) { return getTumorRnaSample(meta).containsKey(FileType.BAM) }


    // Status
    public static boolean hasTumorDna(Map meta) { return hasTumorDnaBam(meta) || hasTumorDnaReduxBam(meta) || hasTumorDnaFastq(meta) }
    public static boolean hasNormalDna(Map meta) { return hasNormalDnaBam(meta) || hasNormalDnaReduxBam(meta) || hasNormalDnaFastq(meta) }
    public static boolean hasDonorDna(Map meta) { return hasDonorDnaBam(meta) || hasDonorDnaReduxBam(meta) || hasDonorDnaFastq(meta) }
    public static boolean hasTumorRna(Map meta) { return hasTumorRnaBam(meta) || hasTumorRnaFastq(meta) }
}
