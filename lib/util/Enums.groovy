package util

class Enums {

    public static String createBulletedList(List<String> items) {
        def bullets = items.collect{ item -> " - ${item}" }
        return bullets.join('\n')
    }

    public static <T extends Enum<T>> String createBulletedList(Class<T> enum_class) {
        def items = enum_class.getEnumConstants().collect{ enumValue -> enumValue.name() }
        return createBulletedList(items)
    }

    public static <T extends Enum<T>> T getValidatedEnumFromString(String string, Class<T> enum_class, boolean ignore_case = true) {

        try {

            def search_string = ignore_case
                ? string.toUpperCase()
                : string

            return Enum.valueOf(enum_class, search_string)

        } catch (IllegalArgumentException e) {

            def enum_class_name = enum_class.getSimpleName()
            def case_sensitive_string = ignore_case ? "(not case sensitive)" : "(case sensitive)"
            def enum_values_string = createBulletedList(enum_class)

            def error_message = "Invalid ${enum_class_name} '${string}'.\n\nValid options ${case_sensitive_string}:\n${enum_values_string}"

            throw new IllegalArgumentException(error_message)
        }
    }

    public static <T extends Enum<T>> T validateEnumFromString(String string, Class<T> enum_class, boolean ignore_case = true) {
        // NOTE(LN): alias method for code clarity
        return getValidatedEnumFromString(string, enum_class, ignore_case)
    }
}