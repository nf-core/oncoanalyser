package pipeline

import util.Enums

public enum PipelineMode {
    PANEL_RESOURCE_CREATION,
    PREPARE_REFERENCE,
    PURITY_ESTIMATE,
    TARGETED,
    WGTS;

    public static PipelineMode fromString(String string) {
        return Enums.getValidatedEnumFromString(string, PipelineMode)
    }
}
