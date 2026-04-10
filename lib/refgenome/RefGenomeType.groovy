package refgenome

public enum RefGenomeType {
    ALT,
    NO_ALT;

    public String getName() { return this.toString().toLowerCase() }

    public static List<String> getNames() {
        return values().collect { type -> type.toString() }
    }
}
