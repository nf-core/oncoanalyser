//
// Process-selection helpers for the nf-core/oncoanalyser pipeline
//

include { Process } from './types'

def getRunStages(processes_include, exclude, manual_select, log) {

    def processes

    if (manual_select) {
        processes = getProcessList(manual_select, log)

        if (processes_include || exclude) {
            log.warn 'When manually selecting processes, including/excluding processes is ignored'
        }

    } else {

        // Get default processes
        processes = Process.values().toList()

        // NOTE(LN): Disable some processes from running by default
        Constants.DEFAULT_EXCLUDED_PROCESSES.each {it -> processes.remove(it) }

        def include_list = getProcessList(processes_include, log)
        def exclude_list = getProcessList(exclude, log)
        checkIncludeExcludeList(include_list, exclude_list, log)

        processes.addAll(include_list)
        processes.removeAll(exclude_list)
    }

    return Process
        .values()
        .collectEntries { p -> [p.name().toLowerCase(), p in processes] }
}

def getProcessList(process_str, log) {
    if (! process_str) {
        return []
    }
    return process_str.tokenize(',').collect { name ->
            try {
                return Process.valueOf(name.toUpperCase())
            } catch(java.lang.IllegalArgumentException e) {
                def processes_str = getProcessNames().join('\n  - ')
                log.error "received invalid process: '${name}'. Valid options are:\n  - ${processes_str}"
                exit 1
            }
        }
        .unique()
}

def checkIncludeExcludeList(include_list, exclude_list, log) {
    def processes_shared = include_list + exclude_list
        .countBy { it }
        .findAll { k, v -> v > 1 }
        .keySet()

    if (processes_shared) {
        def processes_shared_str = processes_shared.join('\n  - ')
        def message_base = 'the following processes were found in the include and the exclude list'
        log.error "${message_base}:\n  - ${processes_shared_str}"
        exit 1
    }
}

def getProcessNames() {
    Process
        .values()
        *.name()
        *.toLowerCase()
}
