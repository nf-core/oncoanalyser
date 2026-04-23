package refdata

import pipeline.PipelineMode
import pipeline.SupportedPanel
import util.Enums
import sample.Inputs

class PrepareReferenceConfig {

    public static Map<String, Boolean> forPipelineRun(List<Map> inputs, PipelineMode pipeline_mode, Map<String, Boolean> stages) {

        def has_dna = inputs.any { Inputs.hasTumorDna(it) }
        def has_rna_fastq = inputs.any { Inputs.hasTumorRnaFastq(it) }
        def has_dna_fastq = inputs.any { Inputs.hasTumorDnaFastq(it) || Inputs.hasNormalDnaFastq(it) }

        return [
            require_fasta: true,
            require_fai: true,
            require_dict: true,
            require_img: true,

            require_bwamem2_index: has_dna_fastq && stages.alignment,
            require_star_index: has_rna_fastq && stages.alignment,

            require_gridss_index: has_dna && pipeline_mode == PipelineMode.WGTS && stages.virusinterpreter,
            require_hmftools_data: true,
            require_panel_data: pipeline_mode == PipelineMode.TARGETED
        ]
    }

    public static Map<String, Boolean> forPrepRefOnly(Map params) {

        def ref_data_types = params.ref_data_types
            .tokenize(',')
            .collect { Enums.getValidatedEnumFromString(it, RefDataType) }

        if (
            ref_data_types.contains(RefDataType.WGS) ||
            ref_data_types.contains(RefDataType.WTS) ||
            ref_data_types.contains(RefDataType.TARGETED)
        ) {
            ref_data_types += [
                RefDataType.FASTA,
                RefDataType.FAI,
                RefDataType.DICT,
                RefDataType.IMG,
                RefDataType.HMFTOOLS
            ]
        }

        if (ref_data_types.contains(RefDataType.WGS)) {
            ref_data_types += [RefDataType.GRIDSS_INDEX]
        }

        if (ref_data_types.contains(RefDataType.TARGETED)) {
            ref_data_types += [RefDataType.PANEL]
        }

        def require_fasta = ref_data_types.contains(RefDataType.FASTA)
        def require_fai = ref_data_types.contains(RefDataType.FAI)
        def require_dict = ref_data_types.contains(RefDataType.DICT)
        def require_img = ref_data_types.contains(RefDataType.IMG)

        def require_bwamem2_index = ref_data_types.contains(RefDataType.BWAMEM2_INDEX) || ref_data_types.contains(RefDataType.DNA_ALIGNMENT)
        def require_star_index = ref_data_types.contains(RefDataType.STAR_INDEX) || ref_data_types.contains(RefDataType.RNA_ALIGNMENT)

        def require_gridss_index = ref_data_types.contains(RefDataType.GRIDSS_INDEX)
        def require_hmftools_data = ref_data_types.contains(RefDataType.HMFTOOLS)
        def require_panel_data = ref_data_types.contains(RefDataType.PANEL)

        if (require_panel_data) {

            if (params.panel == null) {
                throw new IllegalStateException("Preparing panel specific reference data requires the --panel CLI argument to be provided")
            }

            def maybe_supported_panel = SupportedPanel.fromString((String) params.panel)
            if (!maybe_supported_panel) {
                throw new IllegalStateException("Preparing panel specific reference data not supported for custom panel: ${params.panel}")
            }
        }

        return [
            require_fasta: require_fasta,
            require_fai: require_fai,
            require_dict: require_dict,
            require_img: require_img,

            require_bwamem2_index: require_bwamem2_index,
            require_star_index: require_star_index,

            require_gridss_index: require_gridss_index,
            require_hmftools_data: require_hmftools_data,
            require_panel_data: require_panel_data,
        ]
    }

}
