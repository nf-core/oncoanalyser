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

def getSequencingPlatformPons(hmf_data, sequencing_platform_string, log) {
    def sequencing_platform = getEnumFromString(sequencing_platform_string, SequencingPlatform)
    hmf_data.map { d ->
        if (sequencing_platform == SequencingPlatform.ILLUMINA) {
            return [
                'esvee_breakends': d.esvee_pon_breakends_illumina,
                'esvee_breakpoints': d.esvee_pon_breakpoints_illumina,
                'sage': d.sage_pon_illumina,
            ]
        } else if (sequencing_platform == SequencingPlatform.SBX) {
            return [
                'esvee_breakends': d.esvee_pon_breakends_sbx,
                'esvee_breakpoints': d.esvee_pon_breakpoints_sbx,
                'sage': d.sage_pon_sbx,
            ]
        } else if (sequencing_platform == SequencingPlatform.ULTIMA) {
            return [
                'esvee_breakends': d.esvee_pon_breakends_ultima,
                'esvee_breakpoints': d.esvee_pon_breakpoints_ultima,
                'sage': d.sage_pon_ultima,
            ]
        } else {
            log.error "Got bad sequencing platform: ${sequencing_platform}"
            exit 1
        }
    }
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
