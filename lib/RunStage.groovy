import util.Enums

public enum RunStage {

    ALIGNMENT,
    AMBER,
    BAMTOOLS,
    CHORD,
    CIDER,
    COBALT,
    CUPPA,
    ESVEE,
    ISOFOX,
    LILAC,
    LINX,
    NEO,
    ORANGE,
    PAVE,
    PEACH,
    PURPLE,
    QSEE,
    REDUX,
    SAGE,
    SAGE_APPEND,
    SAGE_VIS,
    SIGS,
    TEAL,
    VIRUSINTERPRETER,
    WISP;

    private static List<RunStage> DEFAULT_EXCLUDED_PROCESSES = [] // For experimental tools

    public static Map<String, Boolean> getValidatedRunStages(String include, String exclude, String manual_select, log) {

        def processes

        if (manual_select) {
            processes = getProcessList(manual_select)

            if (include || exclude) {
                log.warning "When manually selecting processes, including/excluding processes is ignored"
            }

        } else {

            // Get default processes
            processes = values().toList()

            // NOTE(LN): Disable some processes from running by default
            DEFAULT_EXCLUDED_PROCESSES.each {it -> processes.remove(it) }

            def include_list = getProcessList(include)
            def exclude_list = getProcessList(exclude)
            checkIncludeExcludeList(include_list, exclude_list, log)

            processes.addAll(include_list)
            processes.removeAll(exclude_list)
        }

        return values().collectEntries { p -> [p.name().toLowerCase(), p in processes] }
    }

    private static List<RunStage> getProcessList(String process_str) {
        if (!process_str) {
            return []
        }
        return process_str
            .tokenize(',')
            .collect { return Enums.getValidatedEnumFromString(it, RunStage) }
            .unique()
    }

    private static void checkIncludeExcludeList(List<RunStage> include_list, List<RunStage> exclude_list, log) {
        def processes_shared = [*include_list, *exclude_list]
            .countBy { it }
            .findAll { k, v -> v > 1 }
            .keySet()

        if (processes_shared) {
            def message_base = 'the following processes was found in the include and the exclude list'
            def processes_shared_str = processes_shared.join('\n  - ')
            throw new IllegalArgumentException("${message_base}:\n  - ${processes_shared_str}")
        }
    }
}
