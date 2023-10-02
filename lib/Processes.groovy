import Constants
import Utils


class Processes {

    public static setProcesses(run_mode, manual_select, targeted_mode, log) {
        def processes = []

        if (manual_select) {
            return processes
        }

        switch(run_mode) {

            case Constants.RunMode.DNA:
                processes = [
                    Constants.Process.AMBER,
                    Constants.Process.BAMTOOLS,
                    Constants.Process.CHORD,
                    Constants.Process.COBALT,
                    Constants.Process.CUPPA,
                    Constants.Process.FLAGSTAT,
                    Constants.Process.GRIDSS,
                    Constants.Process.GRIPSS,
                    Constants.Process.LILAC,
                    Constants.Process.LINX,
                    Constants.Process.ORANGE,
                    Constants.Process.PAVE,
                    Constants.Process.PURPLE,
                    Constants.Process.SAGE,
                    Constants.Process.SIGS,
                    Constants.Process.VIRUSINTERPRETER,
                ]
                break

            case Constants.RunMode.RNA:
                processes = [
                    Constants.Process.CUPPA,
                    Constants.Process.ISOFOX,
                ]
                break

            case Constants.RunMode.DNA_RNA:
                processes = [
                    Constants.Process.AMBER,
                    Constants.Process.BAMTOOLS,
                    Constants.Process.CHORD,
                    Constants.Process.COBALT,
                    Constants.Process.CUPPA,
                    Constants.Process.FLAGSTAT,
                    Constants.Process.GRIDSS,
                    Constants.Process.GRIPSS,
                    Constants.Process.ISOFOX,
                    Constants.Process.LILAC,
                    Constants.Process.LINX,
                    Constants.Process.ORANGE,
                    Constants.Process.PAVE,
                    Constants.Process.PURPLE,
                    Constants.Process.SAGE,
                    Constants.Process.SIGS,
                    Constants.Process.VIRUSINTERPRETER,
                ]
                break

            //case Constants.RunMode.PANEL:
            //    processes = [
            //        Constants.Process.AMBER,
            //        Constants.Process.BAMTOOLS,
            //        Constants.Process.COBALT,
            //        Constants.Process.FLAGSTAT,
            //        Constants.Process.GRIDSS,
            //        Constants.Process.GRIPSS,
            //        Constants.Process.LILAC,
            //        Constants.Process.LINX,
            //        Constants.Process.ORANGE,
            //        Constants.Process.PAVE,
            //        Constants.Process.PURPLE,
            //        Constants.Process.SAGE,
            //    ]
            //    break

            default:
                log.error "\nERROR: we should never have come here"
                System.exit(1)
        }

        return processes
    }

    public static getRunStages(run_mode, include, exclude, manual_select, targeted_mode, log) {
        def processes = this.setProcesses(run_mode, manual_select, targeted_mode, log)
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
                    log.error "\nERROR: recieved invalid process: '${name}'. Valid options are:\n  - ${processes_str}"
                    System.exit(1)
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
            log.error "\nERROR: ${message_base}:\n  - ${processes_shared_str}"
            System.exit(1)
        }
    }

    public static getProcessNames() {
        Constants.Process
            .values()
            *.name()
            *.toLowerCase()
    }
}
