import nextflow.Nextflow

class Enums {

    public static getEnumNames(enum_class, ignore_case = true) {

        def strings = enum_class.values()*.name()

        if(ignore_case) {
            strings = strings*.toLowerCase()
        }

        return strings
    }

    public static getEnumFromString(string, enum_class, ignore_case = true) {
        try {
            def searchString = ignore_case ? string.toUpperCase() : string
            return enum_class.valueOf(searchString)
        } catch(IllegalArgumentException e) {
            return null
        }
    }

    public static getValidatedEnumFromString(string, enum_class, log, ignore_case = true) {

        def enum_value = getEnumFromString(string, enum_class, ignore_case)

        if(enum_value !== null) {
            return enum_value
        } else {
            def enum_class_name = enum_class.getSimpleName()
            def enum_values_string = getEnumNames(enum_class, ignore_case).join('\n  - ')

            log.error "Invalid ${enum_class_name}: '${string}'\n\nValid options are:\n  - ${enum_values_string}"
            Nextflow.exit(1)
        }
    }

    public static validateEnumFromString(string, enum_class, log, ignore_case = true){
        // NOTE(LN): alias method for code clarity
        getValidatedEnumFromString(string, enum_class, log, ignore_case)
    }

}
