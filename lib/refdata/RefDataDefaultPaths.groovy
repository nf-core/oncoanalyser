package refdata

import refgenome.RefGenomeVersion

class RefDataDefaultPaths {

    public static String hmfData(RefGenomeVersion version) {

        return switch (version) {
            case RefGenomeVersion.V37 -> 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.3.0--2.tar.gz'
            case RefGenomeVersion.V38 -> 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.3.0--2.tar.gz'
            default -> throw new IllegalArgumentException("Invalid ref genome version: ${version}")
        }
    }

}
