package panel

public enum Panel {
    TSO500,
    MSK;

    public static Panel fromString(String string) {
        try {
            return valueOf(string)
        } catch (IllegalArgumentException e){
            return null
        }
    }
}
