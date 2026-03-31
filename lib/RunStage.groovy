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

    // Direct dependencies for each stage. The resolver follows the chain transitively,
    // e.g. PURPLE -> PAVE -> SAGE -> REDUX -> ALIGNMENT
    private static Map<RunStage, Set<RunStage>> getDependencies() {
        def S = RunStage
        return [
            (S.ALIGNMENT):         [] as Set,
            (S.REDUX):             [S.ALIGNMENT] as Set,
            (S.BAMTOOLS):          [S.REDUX] as Set,
            (S.ISOFOX):            [S.ALIGNMENT] as Set,
            (S.AMBER):             [S.REDUX] as Set,
            (S.COBALT):            [S.REDUX] as Set,
            (S.ESVEE):             [S.REDUX] as Set,
            (S.SAGE):              [S.REDUX] as Set,
            (S.PAVE):              [S.SAGE] as Set,
            (S.PURPLE):            [S.AMBER, S.COBALT, S.ESVEE, S.PAVE] as Set,
            (S.QSEE):              [S.REDUX, S.BAMTOOLS, S.COBALT, S.ESVEE, S.PURPLE] as Set,
            (S.SAGE_APPEND):       [S.PURPLE] as Set,
            (S.SAGE_VIS):          [S.REDUX, S.PURPLE] as Set,
            (S.LINX):              [S.PURPLE] as Set,
            (S.CHORD):             [S.PURPLE] as Set,
            (S.SIGS):              [S.PURPLE] as Set,
            (S.PEACH):             [S.PURPLE] as Set,
            (S.LILAC):             [S.REDUX, S.PURPLE] as Set,
            (S.TEAL):              [S.REDUX, S.BAMTOOLS, S.COBALT, S.PURPLE] as Set,
            (S.VIRUSINTERPRETER):  [S.REDUX, S.BAMTOOLS, S.PURPLE] as Set,
            (S.CIDER):             [S.REDUX] as Set,
            (S.CUPPA):             [S.PURPLE, S.LINX, S.ISOFOX, S.VIRUSINTERPRETER] as Set,
            (S.NEO):               [S.PURPLE, S.ISOFOX, S.LILAC, S.LINX, S.SAGE_APPEND] as Set,
            (S.ORANGE):            [S.SAGE, S.SAGE_APPEND, S.PURPLE, S.QSEE, S.LINX, S.VIRUSINTERPRETER, S.CHORD, S.SIGS, S.LILAC, S.CUPPA, S.PEACH, S.ISOFOX] as Set,
            (S.WISP):              [S.AMBER, S.COBALT, S.SAGE_APPEND] as Set,
        ]
    }

    // Maps samplesheet file types to the stage(s) they satisfy.
    // A BAM satisfies ALIGNMENT; a BAM_REDUX satisfies both ALIGNMENT and REDUX.
    private static Map<SampleMeta.FileType, Set<RunStage>> getInputSatisfiesStage() {
        def S = RunStage
        def F = SampleMeta.FileType
        return [
            (F.BAM):                  [S.ALIGNMENT] as Set,
            (F.CRAM):                 [S.ALIGNMENT] as Set,
            (F.BAM_REDUX):            [S.ALIGNMENT, S.REDUX] as Set,
            (F.CRAM_REDUX):           [S.ALIGNMENT, S.REDUX] as Set,
            (F.AMBER_DIR):            [S.AMBER] as Set,
            (F.COBALT_DIR):           [S.COBALT] as Set,
            (F.ESVEE_DIR):            [S.ESVEE] as Set,
            (F.SAGE_DIR):             [S.SAGE] as Set,
            (F.PAVE_DIR):             [S.PAVE] as Set,
            (F.PURPLE_DIR):           [S.PURPLE] as Set,
            (F.ISOFOX_DIR):           [S.ISOFOX] as Set,
            (F.BAMTOOLS_DIR):         [S.BAMTOOLS] as Set,
            (F.LINX_ANNO_DIR):        [S.LINX] as Set,
            (F.LILAC_DIR):            [S.LILAC] as Set,
            (F.CHORD_DIR):            [S.CHORD] as Set,
            (F.CUPPA_DIR):            [S.CUPPA] as Set,
            (F.PEACH_DIR):            [S.PEACH] as Set,
            (F.QSEE_DIR):             [S.QSEE] as Set,
            (F.SAGE_APPEND_DIR):      [S.SAGE_APPEND] as Set,
            (F.SIGS_DIR):             [S.SIGS] as Set,
            (F.VIRUSINTERPRETER_DIR): [S.VIRUSINTERPRETER] as Set,
        ]
    }

    public static Map<String, Boolean> getValidatedRunStages(String include, String exclude, String manual_select, List inputs, log) {

        def processes

        if (manual_select) {
            processes = getProcessList(manual_select)

            if (include || exclude) {
                log.warning "When manually selecting processes, including/excluding processes is ignored"
            }

        } else if (include) {

            def target_list = getProcessList(include)
            def exclude_list = getProcessList(exclude)
            def satisfied_stages = getSatisfiedStages(inputs)

            def resolved = resolveDependencies(target_list as Set, satisfied_stages)
            validateExclusions(resolved, exclude_list, target_list)
            resolved.removeAll(exclude_list)

            logResolvedStages(target_list, resolved, satisfied_stages, exclude_list, log)

            processes = resolved.toList()

        } else {

            // Default: all stages enabled minus exclude list
            processes = values().toList()

            // NOTE(LN): Disable some processes from running by default
            DEFAULT_EXCLUDED_PROCESSES.each { it -> processes.remove(it) }

            def exclude_list = getProcessList(exclude)
            processes.removeAll(exclude_list)
        }

        return values().collectEntries { p -> [p.name().toLowerCase(), p in processes] }
    }

    private static Set<RunStage> resolveDependencies(Set<RunStage> targets, Set<RunStage> satisfied_stages) {
        def required = [] as Set
        def queue = new ArrayDeque<RunStage>(targets)

        while (!queue.isEmpty()) {
            def stage = queue.poll()

            if (stage in required) {
                continue
            }

            required.add(stage)

            // If this stage is satisfied by user input, don't traverse its upstream dependencies
            if (stage in satisfied_stages) {
                continue
            }

            def deps = getDependencies()[stage] ?: ([] as Set)
            for (dep in deps) {
                if (!(dep in required)) {
                    queue.add(dep)
                }
            }
        }

        // Remove stages satisfied by input — they don't need to run
        return required - satisfied_stages
    }

    private static Set<RunStage> getSatisfiedStages(List inputs) {
        if (!inputs) {
            return [] as Set
        }

        def input_satisfies = getInputSatisfiesStage()

        // Collect satisfied stages per group, then intersect across all groups.
        // A stage is only globally satisfied if ALL groups provide the pre-computed input.
        def per_group_satisfied = inputs.collect { meta ->
            def group_satisfied = [] as Set
            SampleMeta.INPUT.each { input_key, input_def ->
                def (file_type, sample_types, sequence_type) = input_def
                if (input_satisfies.containsKey(file_type)) {
                    if (Inputs.hasExistingInput(meta, input_def)) {
                        group_satisfied.addAll(input_satisfies[file_type])
                    }
                }
            }
            return group_satisfied
        }

        // Intersect: only keep stages satisfied in ALL groups
        def result = per_group_satisfied[0] as Set
        per_group_satisfied.each { group_satisfied ->
            result.retainAll(group_satisfied)
        }

        return result
    }

    private static void validateExclusions(Set<RunStage> resolved_stages, List<RunStage> exclude_list, List<RunStage> targets) {
        def target_set = targets as Set
        def exclude_set = exclude_list as Set

        // Check if user is trying to exclude a stage they explicitly requested
        def excluded_targets = target_set.intersect(exclude_set)
        if (excluded_targets) {
            def excluded_targets_str = excluded_targets.collect { it.name().toLowerCase() }.join(', ')
            throw new IllegalArgumentException(
                "Cannot both include and exclude the same stage: [${excluded_targets_str}]."
            )
        }

        // Check if user is trying to exclude a required dependency
        def conflicts = resolved_stages.intersect(exclude_set)
        if (conflicts) {
            def conflicts_str = conflicts.collect { it.name().toLowerCase() }.join(', ')
            def targets_str = targets.collect { it.name().toLowerCase() }.join(', ')
            throw new IllegalArgumentException(
                "Cannot exclude [${conflicts_str}] because required as dependencies of [${targets_str}].\n" +
                "To skip these stages, provide their pre-computed outputs in the samplesheet instead."
            )
        }
    }

    private static void logResolvedStages(List<RunStage> targets, Set<RunStage> resolved, Set<RunStage> satisfied, List<RunStage> excluded, log) {
        def targets_str = targets.collect { it.name().toLowerCase() }.join(', ')
        def resolved_str = resolved.collect { it.name().toLowerCase() }.sort().join(', ')
        def satisfied_str = satisfied ? satisfied.collect { it.name().toLowerCase() }.sort().join(', ') : '-'
        def excluded_str = excluded ? excluded.collect { it.name().toLowerCase() }.sort().join(', ') : '-'

        log.info "Dependency resolution for processes_include [${targets_str}]:"
        log.info "  Stages to run:            ${resolved_str}"
        log.info "  Skipped (input provided):  ${satisfied_str}"
        log.info "  Skipped (excluded):        ${excluded_str}"
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
}
