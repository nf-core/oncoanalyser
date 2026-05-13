package util

class Messages {

    public static void error(String... message){

        def message_string = "\n"
        message_string += message.join('\n')

        throw new RuntimeException(message_string)
    }

    public static String createBulletedList(List<String> items) {
        def bullets = items.collect{ item -> " - ${item}" }
        return bullets.join('\n')
    }

    public static <T extends Enum<T>> String createBulletedList(Class<T> enum_class) {
        def items = enum_class.getEnumConstants().collect{ enumValue -> enumValue.name() }
        return Messages.createBulletedList(items)
    }
}