import nextflow.Nextflow

import Constants
import Utils


class Processes {

    public static getRunStages(include, exclude, manual_select, log) {

        // Get default processes
        // NOTE(SW): currently set all except Neo to run by default; Process.NEO excluded to be more concise in code
        def processes
        if (manual_select) {
            processes = []
        } else {
            processes = Constants.Process.values().toList()
            processes.remove(Constants.Process.NEO)
        }

        def include_list = this.getProcessList(include, log)
        def exclude_list = this.getProcessList(exclude, log)
        this.checkIncludeExcludeList(include_list, exclude_list, log)

        processes.addAll(include_list)
        processes.removeAll(exclude_list)

        return Constants.Process
            .values()
            .collectEntries { p -> [p.name().toLowerCase(), p in processes] }
    }

    public static getProcessList(process_str, log) {
        if (!process_str) {
            return []
        }
        return process_str
            .tokenize(',')
            .collect { name ->
                try {
                    return Constants.Process.valueOf(name.toUpperCase())
                } catch(java.lang.IllegalArgumentException e) {
                    def processes_str = Processes.getProcessNames().join('\n  - ')
                    log.error "received invalid process: '${name}'. Valid options are:\n  - ${processes_str}"
                    Nextflow.exit(1)
                }
            }
            .unique()
    }

    public static checkIncludeExcludeList(include_list, exclude_list, log) {
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

    public static getProcessNames() {
        Constants.Process
            .values()
            *.name()
            *.toLowerCase()
    }
}
