import nextflow.Nextflow

import Constants
import Utils


class Processes {

    public static getValidatedRunStages(include, exclude, manual_select, log) {

        def processes

        if (manual_select) {
            processes = getProcessList(manual_select, log)

            if (include || exclude) {
                log.warning "When manually selecting processes, including/excluding processes is ignored"
            }

        } else {

            // Get default processes
            processes = Constants.Process.values().toList()

            // NOTE(LN): Disable some processes from running by default
            Constants.DEFAULT_EXCLUDED_PROCESSES.each {it -> processes.remove(it) }

            def include_list = getProcessList(include, log)
            def exclude_list = getProcessList(exclude, log)
            checkIncludeExcludeList(include_list, exclude_list, log)

            processes.addAll(include_list)
            processes.removeAll(exclude_list)
        }

        return Constants.Process
            .values()
            .collectEntries { p -> [p.name().toLowerCase(), p in processes] }
    }

    private static getProcessList(process_str, log) {
        if (!process_str) {
            return []
        }
        return process_str
            .tokenize(',')
            .collect { return Utils.getValidatedEnumFromString(it, Constants.Process, log) }
            .unique()
    }

    private static checkIncludeExcludeList(include_list, exclude_list, log) {
        def processes_shared = [*include_list, *exclude_list]
            .countBy { it }
            .findAll { k, v -> v > 1 }
            .keySet()

        if (processes_shared) {
            def processes_shared_str = processes_shared.join('\n  - ')
            def message_base = 'the following processes was found in the include and the exclude list'
            log.error "${message_base}:\n  - ${processes_shared_str}"
            Nextflow.exit(1)
        }
    }
}
