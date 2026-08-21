//
// Generic helpers for the nf-core/oncoanalyser pipeline
//

include { RunMode            } from './types'
include { SequencingPlatform } from './types'

def createStubPlaceholders(params) {

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

    if (params.panel != null) {
        params.panel_data_paths[params.panel.toLowerCase()][params.genome_version.toString()]
            .each { k, v ->
                fps << "${params.ref_data_panel_data_path.replaceAll('/$', '')}/${v}"
            }
    }

    fps.each { fp_str ->
        if (fp_str == null) {
            return
        }

        def fp = getFileObject(fp_str)

        if (! fp_str || fp.exists()) {
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

def getSequencingPlatformPon(hmf_data, sequencing_platform_string, pon_key, log) {
    def sequencing_platform = getEnumFromString(sequencing_platform_string, SequencingPlatform)
    def suffix = null
    if (sequencing_platform == SequencingPlatform.ILLUMINA) {
        suffix = 'illumina'
    } else if (sequencing_platform == SequencingPlatform.SBX) {
        suffix = 'sbx'
    } else if (sequencing_platform == SequencingPlatform.ULTIMA) {
        suffix = 'ultima'
    } else {
        log.error "Got bad sequencing platform: ${sequencing_platform}"
        exit 1
    }

    def pon_property = pon_key == 'esvee_breakends' ? "esvee_pon_breakends_${suffix}"
        : pon_key == 'esvee_breakpoints' ? "esvee_pon_breakpoints_${suffix}"
        : "sage_pon_${suffix}"

    // NOTE(SW): map over hmf_data (a broadcast) rather than a shared intermediate, so each
    // call yields an independent stream a typed workflow may consume once.
    hmf_data.map { d -> d[pon_property] }
}

def getEnumFromString(s, e) {
    try {
        return e.valueOf(s.toUpperCase())
    } catch (java.lang.IllegalArgumentException err) {
        return null
    }
}

def getEnumNames(e) {
    e
        .values()
        *.name()
        *.toLowerCase()
}

def getEnumFromStringOrFail(value, enum_class, label, log) {
    def enum_value = getEnumFromString(value, enum_class)
    if (! enum_value) {
        def options = getEnumNames(enum_class).join('\n  - ')
        log.error "received invalid ${label}: '${value}'. Valid options are:\n  - ${options}"
        exit 1
    }
    return enum_value
}

def getFileObject(path) {
    return path ? nextflow.Nextflow.file(path) : null
}

def getRunMode(run_mode, log) {
    def run_mode_enum = getEnumFromString(run_mode, RunMode)
    if (! run_mode_enum) {
        def run_modes_str = getEnumNames(RunMode).join('\n  - ')
        log.error "received an invalid run mode: '${run_mode}'. Valid options are:\n  - ${run_modes_str}"
        exit 1
    }
    return run_mode_enum
}

def selectCurrentOrExisting(val, existing) {
    return existing != null ? existing : val
}
