//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import nextflow.Nextflow

class Utils {

    public static void createStubPlaceholders(params) {

        def fps = [
            params.ref_data_genome_alt,
            params.ref_data_genome_bwamem2_index,
            params.ref_data_genome_dict,
            params.ref_data_genome_fai,
            params.ref_data_genome_fasta,
            params.ref_data_genome_gridss_index,
            params.ref_data_genome_gtf,
            params.ref_data_genome_star_index,
        ]

        params.hmf_data_paths[params.genome_version.toString()]
            .each { k, v ->
                fps << "${params.ref_data_hmf_data_path.replaceAll('/$', '')}/${v}"
            }

        if (params.panel !== null) {
            params.panel_data_paths[params.panel][params.genome_version.toString()]
                .each { k, v ->
                    fps << "${params.ref_data_panel_data_path.replaceAll('/$', '')}/${v}"
                }
        }

        fps.each { fp_str ->
            if (fp_str === null) {
                return
            }

            def fp = getFileObject(fp_str)

            if (!fp_str || fp.exists()) {
                return
            }

            if (fp_str.endsWith('/')) {
                fp.mkdirs()
            } else {
                fp.getParent().mkdirs()
                fp.toFile().createNewFile()
            }
        }
    }

    public static getFileObject(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    // Enums
    public static getEnumNames(enum_class, ignore_case = true) {

        def strings = enum_class.values()*.name()

        if(ignore_case) {
            strings = strings*.toLowerCase()
        }

        return strings
    }

    public static getEnumFromString(string, enum_class, ignore_case = true) {
        try {
            def searchString = ignore_case ? string.toUpperCase() : string
            return enum_class.valueOf(searchString)
        } catch(IllegalArgumentException e) {
            return null
        }
    }

    public static getValidatedEnumFromString(string, enum_class, log, ignore_case = true) {

        def enum_value = getEnumFromString(string, enum_class, ignore_case)

        if(enum_value !== null) {
            return enum_value
        } else {
            def enum_class_name = enum_class.getSimpleName()
            def enum_values_string = getEnumNames(enum_class, ignore_case).join('\n  - ')

            log.error "Invalid ${enum_class_name}: '${string}'\n\nValid options are:\n  - ${enum_values_string}"
            Nextflow.exit(1)
        }
    }

    public static validateEnumFromString(string, enum_class, log, ignore_case = true){
        // NOTE(LN): alias method for code clarity
        getValidatedEnumFromString(string, enum_class, log, ignore_case)
    }

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

    public static getTumorDnaBai(meta) {
        return getTumorDnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
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
    public static getNormalDnaBai(meta) {
        return getNormalDnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
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

    public static hasDnaFastq(meta) {
        return hasNormalDnaFastq(meta) || hasTumorDnaFastq(meta)
    }

    public static hasDnaReduxBam(meta) {
        return hasNormalDnaReduxBam(meta) || hasTumorDnaReduxBam(meta)
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

    public static getDonorDnaBai(meta) {
        return getDonorDnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
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

    public static getTumorRnaBai(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
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
