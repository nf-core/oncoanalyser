# nf-core/oncoanalyser: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project mostly adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[2.2.0](https://github.com/nf-core/oncoanalyser/releases/tag/2.2.0)] Royal Spoonbill - 2025-08-02

- [241](https://github.com/nf-core/oncoanalyser/pull/241) - Apply minor fixes and updates
  - Allow 'prepare reference' feature to be driven by samplesheet
  - Set minimum stringr / stringi version for CUPPA environment
  - Reintroduce decoy sequences for ESVEE with GRCh37 genomes
  - Update WiGiTS reference data paths
  - Improve FASTQ and longitudinal sample input handling
  - Fix REDUX TSV collection in SAGE append subworkflow
  - Update CHANGELOG.md
- [235](https://github.com/nf-core/oncoanalyser/pull/235) - Publish selected command / log files
- [234](https://github.com/nf-core/oncoanalyser/pull/234) - Update WiGiTS tools and reference data
- [233](https://github.com/nf-core/oncoanalyser/pull/233) - Update documentation
- [232](https://github.com/nf-core/oncoanalyser/pull/232) - Extend the 'prepare reference' functionality
- [231](https://github.com/nf-core/oncoanalyser/pull/231) - Implement 'purity estimate' (WISP) workflow
- [230](https://github.com/nf-core/oncoanalyser/pull/230) - Implement 'panel resource creation' workflow
- [220](https://github.com/nf-core/oncoanalyser/pull/220) - Add reports to tower.yml
- [222](https://github.com/nf-core/oncoanalyser/pull/222) - Post-release bump

### Software dependencies

| Dependency         | Old version | New version |
| ------------------ | ----------- | ----------- |
| `AMBER`            | 4.1.1       | 4.2         |
| `BamTools`         | 1.3         | 1.4.2       |
| `bwa-plus`         | 1.0.0       | -           |
| `bwa-mem2`         | -           | 2.3         |
| `CHORD`            | 2.1.0       | 2.1.2       |
| `COBALT`           | 2.0         | 2.1         |
| `ESVEE`            | 1.0.3       | 1.1.2       |
| `ISOFOX`           | 1.7.1       | 1.7.2       |
| `LILAC`            | 1.6         | 1.7.1       |
| `LINX`             | 2.0.2       | 2.1         |
| `NEO`              | 1.2         | 1.2.1       |
| `ORANGE`           | 3.8.1       | 4.1         |
| `PAVE`             | 1.7.1       | 1.8         |
| `PURPLE`           | 4.1         | 4.2         |
| `REDUX`            | 1.1.2       | 1.2         |
| `SAGE`             | 4.0         | 4.1         |
| `VirusInterpreter` | 1.7         | 1.7.1       |
| `WISP`             | -           | 1.2         |

### Reference data

| Name                     | Old version | New version |
| ------------------------ | ----------- | ----------- |
| `HMF pipeline resources` | `2.1.0--1`  | `2.2.0--3`  |
| `HMF TSO500 resources`   | `2.0.0--3`  | `2.2.0--3`  |

### Parameters

| Old name    | New name               |
| ----------- | ---------------------- |
| `fastp_umi` | `fastp_umi_enabled`    |
| `redux_umi` | `redux_umi_enabled`    |
| -           | `purity_estimate_mode` |
| -           | `ref_data_types`       |
| -           | `driver_gene_panel`    |
| -           | `target_regions_bed`   |
| -           | `hmftools_log_level`   |

## [[2.1.0](https://github.com/nf-core/oncoanalyser/releases/tag/2.1.0)] Peruvian Pelican - 2025-06-30

- [219](https://github.com/nf-core/oncoanalyser/pull/219) - Add metromap-style diagram for pipeline overview
- [217](https://github.com/nf-core/oncoanalyser/pull/217) - Fix targeted mode parameters
- [210](https://github.com/nf-core/oncoanalyser/pull/210) - Switch to bwa-plus from bwa-mem2
- [207](https://github.com/nf-core/oncoanalyser/pull/207) - Apply minor fixes and updates
  - Update the HMF SAGE PON
  - Bump HMF pipeline resource bundles to 2.1.0--1
  - Bump LINX to 2.0.2
  - Bump PAVE to 1.7.1
  - Update CHANGELOG.md
- [205](https://github.com/nf-core/oncoanalyser/pull/205) - Fix output prepared reference data directory names
- [189](https://github.com/nf-core/oncoanalyser/pull/189) - Implement TEAL subworkflow
- [188](https://github.com/nf-core/oncoanalyser/pull/188) - Implement CIDER subworkflow
- [187](https://github.com/nf-core/oncoanalyser/pull/187) - Implement PEACH subworkflow
- [201](https://github.com/nf-core/oncoanalyser/pull/201) - Fix REDUX TSV discovery for non-local files
- [190](https://github.com/nf-core/oncoanalyser/pull/190) - Fix integer overflow of fastp `--split_by_lines` value
- [199](https://github.com/nf-core/oncoanalyser/pull/199) - Post-release bump

## [[2.0.0](https://github.com/nf-core/oncoanalyser/releases/tag/2.0.0)] Flame Robin - 2025-04-10

- [178](https://github.com/nf-core/oncoanalyser/pull/178) - Apply minor fixes and updates
  - Bump REDUX to 1.1.2
  - Bump ESVEE to 1.0.3
  - Bump ORANGE to 3.8.1
  - Update Hartwig pipeline and TSO500 panel resources
  - Restore missing fastp threads argument
  - Update CHANGELOG.md
- [176](https://github.com/nf-core/oncoanalyser/pull/176) - Fix utilisation of user-provided BamTools directory
- [174](https://github.com/nf-core/oncoanalyser/pull/174) - Fix check for existing LINX inputs
- [173](https://github.com/nf-core/oncoanalyser/pull/173) - Fix BAM / CRAM index discovery on cloud platforms
- [172](https://github.com/nf-core/oncoanalyser/pull/172) - Require Virus Interpreter results to run CUPPA in DNA mode
- [171](https://github.com/nf-core/oncoanalyser/pull/171) - Pin VIRUSBreakend dependency RepeatMasker to 4.1.5
- [170](https://github.com/nf-core/oncoanalyser/pull/170) - Bump linxreport to 1.1.0
- [169](https://github.com/nf-core/oncoanalyser/pull/169) - Bump BWA to 0.7.19
- [168](https://github.com/nf-core/oncoanalyser/pull/168) - Switch to the new Hartwig GRCh38 ALT masked reference genome
- [167](https://github.com/nf-core/oncoanalyser/pull/167) - Disallow CRAM RNA input
- [160](https://github.com/nf-core/oncoanalyser/pull/160) - Fix prepare reference subworkflow
- [159](https://github.com/nf-core/oncoanalyser/pull/159) - Restore CUPPA RNA-only mode
- [158](https://github.com/nf-core/oncoanalyser/pull/158) - Prevent CUPPA writing into the Isofox work directory, invalidating cache
- [75](https://github.com/nf-core/oncoanalyser/pull/75) - Move fastp trimming arguments to modules.config
- [157](https://github.com/nf-core/oncoanalyser/pull/157) - Update ORANGE and Neo entrypoints
- [156](https://github.com/nf-core/oncoanalyser/pull/156) - Apply JVM heap space modifier to new processes
- [155](https://github.com/nf-core/oncoanalyser/pull/155) - Qualify process names in subworkflows to disambigutate
- [141](https://github.com/nf-core/oncoanalyser/pull/141) - Make version collection more robust
- [154](https://github.com/nf-core/oncoanalyser/pull/154) - Update meta.yml, environment.yml and misc generic tools
- [149](https://github.com/nf-core/oncoanalyser/pull/149) - Fix BamTools and ESVEE stubs
- [148](https://github.com/nf-core/oncoanalyser/pull/148) - Implement WiGiTS 2.0.0 (fka WiGiTS 6.0)
- [126](https://github.com/nf-core/oncoanalyser/pull/126) - Expose ORANGE cancer type parameter
- [117](https://github.com/nf-core/oncoanalyser/pull/117) - Add custom panel support
- [116](https://github.com/nf-core/oncoanalyser/pull/116) - Implement Neo subworkflow
- [115](https://github.com/nf-core/oncoanalyser/pull/115) - Add donor sample support
- [107](https://github.com/nf-core/oncoanalyser/pull/107) - Resource requirement tweaks
- [105](https://github.com/nf-core/oncoanalyser/pull/105) - Adjust MarkDups outputs to improve efficiency with k8s
- [98](https://github.com/nf-core/oncoanalyser/pull/98) - Fix typos in error messages for process and run mode check
- [96](https://github.com/nf-core/oncoanalyser/pull/96) - Add missing type field to an entry in the nextflow_schema.json
- [95](https://github.com/nf-core/oncoanalyser/pull/95) - Post-release bump

## [[1.0.0](https://github.com/nf-core/oncoanalyser/releases/tag/1.0.0)] Pied Currawong - 2024-08-26

Initial release of nf-core/oncoanalyser, created with the [nf-core](https://nf-co.re/) template.
