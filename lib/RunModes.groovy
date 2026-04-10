import util.Enums

class RunModes {

    public static final Integer DEFAULT_ISOFOX_READ_LENGTH_WTS = 151
    public static final Integer DEFAULT_ISOFOX_READ_LENGTH_TARGETED = 93

    public static enum Pipeline {
        PANEL_RESOURCE_CREATION,
        PREPARE_REFERENCE,
        PURITY_ESTIMATE,
        TARGETED,
        WGTS;

        public static Pipeline fromString(String string) {
            return Enums.getValidatedEnumFromString(string, Pipeline)
        }
    }

    public static enum PurityEstimate {
        TARGETED,
        WGTS;

        public static PurityEstimate fromString(String string) {
            return Enums.getValidatedEnumFromString(string, PurityEstimate)
        }
    }

    public static enum SequencingType {
        ILLUMINA,
        SBX,
        ULTIMA;

        public static SequencingType fromString(String string) {
            return Enums.getValidatedEnumFromString(string, SequencingType, false)
        }
    }
}
