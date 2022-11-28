import Constants


class Processes {

  public static set_processes(mode_str, log) {
    def mode_enum

    try {
      mode_enum = Constants.PipelineMode.valueOf(mode_str.toUpperCase())
    } catch(java.lang.IllegalArgumentException e) {
      def workflows_str = Processes.get_workflow_names().join('\n  - ')
      log.error "\nERROR: recieved invalid pipeline mode: '${mode_str}'. Valid options are:\n  - ${workflows_str}"
      System.exit(1)
    }

    def processes = []
    switch(mode_enum) {
      case Constants.PipelineMode.FULL:
        processes = Constants.Process.values() as List
        break
      case Constants.PipelineMode.MANUAL:
        break
      case Constants.PipelineMode.GRIDSS_PURPLE_LINX:
        processes = [
          Constants.Process.AMBER,
          Constants.Process.COBALT,
          Constants.Process.GRIDSS,
          Constants.Process.GRIPSS,
          Constants.Process.LINX,
          Constants.Process.PURPLE,
          Constants.Process.SVPREP,
        ]
        break
      case Constants.PipelineMode.CUPPA:
        processes = [
          Constants.Process.AMBER,
          Constants.Process.COBALT,
          Constants.Process.COLLECTWGSMETRICS,
          Constants.Process.CUPPA,
          Constants.Process.GRIDSS,
          Constants.Process.GRIPSS,
          Constants.Process.ISOFOX,
          Constants.Process.LINX,
          Constants.Process.PURPLE,
          Constants.Process.SVPREP,
          Constants.Process.VIRUSINTERPRETER,
        ]
        break
      default:
        log.error "\nERROR: we should never have come here"
        System.exit(1)
    }
    return processes
  }

  public static get_process_list(process_str, log) {
    return process_str
      .tokenize(',')
      .collect { name ->
        try {
          return Constants.Process.valueOf(name.toUpperCase())
        } catch(java.lang.IllegalArgumentException e) {
          def processes_str = Processes.get_process_names().join('\n  - ')
          log.error "\nERROR: recieved invalid process: '${name}'. Valid options are:\n  - ${processes_str}"
          System.exit(1)
        }
      }
      .unique()
  }

  public static check_include_exclude_list(include_list, exclude_list, log) {
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

  public static get_workflow_names() {
    Constants.PipelineMode
      .values()
      *.name()
      *.toLowerCase()
  }

  public static get_process_names() {
    Constants.Process
      .values()
      *.name()
      *.toLowerCase()
  }

}
