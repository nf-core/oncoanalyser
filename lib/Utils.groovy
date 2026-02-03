//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

class Utils {

    // Sample records
    public static getTumorDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.DNA], [:])
    }

    public static getTumorRnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.RNA], [:])
    }

    public static getNormalDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.NORMAL, Constants.SequenceType.DNA], [:])
    }

    public static getDonorDnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.DONOR, Constants.SequenceType.DNA], [:])
    }

    // Sample names
    public static getTumorDnaSampleName(Map named_args, meta) {
        def meta_sample = getTumorDnaSample(meta)
        def sample_id

        if (named_args.getOrDefault('primary', false)) {
            sample_id = meta_sample['sample_id']
        } else {
            sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])
        }

        return sample_id
    }

    public static getTumorDnaSampleName(meta) {
        getTumorDnaSampleName([:], meta)
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
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getTumorDnaBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    public static getTumorDnaReduxBam(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAM_REDUX, null)
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
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getNormalDnaBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    public static getNormalDnaReduxBam(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAM_REDUX, null)
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
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getDonorDnaBam(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    public static getDonorDnaReduxBam(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.BAM_REDUX, null)
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
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    public static getTumorRnaBam(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
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

    // Files - Tool outputs
    //
    // NOTE(LN): We construct the output file paths based on the tool output dir. This so that the user can resume a run by providing
    // for example purple_dir in the sample sheet.
    //
    // Processes in theory could emit these output files, and downstream processes could consume the output file channels. However,
    // this would make resuming messy as multiple e.g. PURPLE output files would need to be provided in the sample sheet

    public static getPurpleSomaticVcf(meta, purple_dir) {
        return nextflow.Nextflow.file(purple_dir).resolve("${getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz")
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

    // Misc
    public static getInput(meta, key) {

        def result = []
        def (key_filetype, key_filetypes, key_sequencetypes) = key

        for (key_sample in [key_filetypes, key_sequencetypes].combinations()) {
            if (meta.containsKey(key_sample) && meta[key_sample].containsKey(key_filetype)) {
                result = meta[key_sample].get(key_filetype)
                break
            }
        }
        return result
    }

    public static hasExistingInput(meta, key) {
        return getInput(meta, key) != []
    }

    public static selectCurrentOrExisting(val, meta, key) {
        if (hasExistingInput(meta, key)) {
            return getInput(meta, key)
        } else {
            return val
        }
    }

}
