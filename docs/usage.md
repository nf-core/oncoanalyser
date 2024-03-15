# nf-core/oncoanalyser: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/oncoanalyser/usage](https://nf-co.re/oncoanalyser/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The oncoanalyser pipeline typically runs from FASTQs or BAMs and supports two modes: (1) whole genome and/or
transcriptome, and (2) targeted panel. Launching an analysis requires only the creation of a samplesheet that describes
details of each input such as the sample type (tumor or normal), sequence type (DNA or RNA), and filepath.

Various aspects of an oncoanalyser analysis can be configured to fit a range of needs, and many of these are considered
[advanced usage](#advanced-usage) of the pipeline. The most useful include:

- precise process selection
- starting from existing data
- granular control over resource data

These features enable oncoanalyser to be run in a highly flexible way. For example, an analysis can be run with existing
PURPLE data as the starting point and skip variant calling processes. Additionally, resource/reference data can staged
locally to optimise execution or modified to create user-defined driver gene panels.

> [!WARNING]
> There are important requirements when using BAMs as input instead of FASTQs:
>
> - STAR must have been run with [specific
>   parameters](https://github.com/hartwigmedical/hmftools/tree/master/isofox#a-note-on-alignment-and-multi-mapping),
>   this is critical for WTS data, and
> - reads are expected to have been aligned to one of the Hartwig-distributed reference genomes (user-defined genomes may be used though are not recommended)

## Supported analyses

A variety of analyses are accessible in oncoanalyser and are implicitly run according to the data available in the
samplesheet. The supported analysis types for each workflow are listed below.

| Input sequence data                 |  WGS/WTS workflow  | Targeted sequencing workflow<sup>\*</sup> |
| ----------------------------------- | :----------------: | :---------------------------------------: |
| • Tumor/normal DNA<br />• Tumor RNA | :white_check_mark: |                     -                     |
| • Tumor only DNA<br />• Tumor RNA   | :white_check_mark: |            :white_check_mark:             |
| • Tumor/normal DNA                  | :white_check_mark: |                     -                     |
| • Tumor only DNA                    | :white_check_mark: |            :white_check_mark:             |
| • Tumor only RNA                    | :white_check_mark: |                     -                     |

<sub><sup>\*</sup> Supported analyses relate to the TSO500 panel only</sub>

## Samplesheet

A samplesheet that contains information of each input in CSV format is needed to run oncoanalyser. The required input
details and columns are [described below](#column-descriptions).

The oncoanalyser pipeline also recognises several input filetypes, including intermediate output files generated during
execution such as the PURPLE output directory. The full list recognised input filetypes is available
[here](https://github.com/nf-core/oncoanalyser/blob/v0.3.1/lib/Constants.groovy#L56-L86).

### Simple example

#### FASTQ

> [!NOTE]
> Currently only non-interleaved paired-end reads are accepted as FASTQ input

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
P1__wgts,P1,SA,normal,dna,fastq,library_id:SA_library;lane:001,/path/to/P1.SA.normal.dna.wgs.001.R1.fastq.gz;/path/to/P1.SA.normal.dna.wgs.001.R2.fastq.gz
P1__wgts,P1,SB,tumor,dna,fastq,library_id:SB_library;lane:001,/path/to/P1.SB.tumor.dna.wgs.001.R1.fastq.gz;/path/to/P1.SB.tumor.dna.wgs.001.R2.fastq.gz
P1__wgts,P1,SC,tumor,rna,fastq,library_id:SC_library;lane:001,/path/to/P1.SC.tumor.rna.wts.001.R1.fastq.gz;/path/to/P1.SC.tumor.rna.wts.001.R2.fastq.gz
```

#### BAM

> [!NOTE]
> Inputs with the `bam` filetype will be processed by MarkDups as required by hmftools. Where an input BAM has already
> been processed by MarkDups, you can avoid needless reprocessing by setting `bam_markdups` as the filetype instead.
>
> Please note there are important requirements around the use of BAMs, see the warning above in the
> [Introduction](#introduction).

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
P1__wgts,P1,SA,normal,dna,bam,/path/to/P1.SA.normal.dna.wgs.bam
P1__wgts,P1,SB,tumor,dna,bam,/path/to/P1.SB.tumor.dna.wgs.bam
P1__wgts,P1,SC,tumor,rna,bam,/path/to/P1.SC.tumor.rna.wts.bam
```

### Multiple lanes

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
P1__wgts,P1,SA,normal,dna,fastq,library_id:SA_library;lane:001,/path/to/P1.SA.normal.dna.wgs.001.R1.fastq.gz;/path/to/P1.SA.normal.dna.wgs.001.R2.fastq.gz
P1__wgts,P1,SA,normal,dna,fastq,library_id:SA_library;lane:002,/path/to/P1.SA.normal.dna.wgs.002.R1.fastq.gz;/path/to/P1.SA.normal.dna.wgs.002.R2.fastq.gz
P1__wgts,P1,SB,tumor,dna,fastq,library_id:SB_library;lane:001,/path/to/P1.SB.tumor.dna.wgs.001.R1.fastq.gz;/path/to/P1.SB.tumor.dna.wgs.001.R2.fastq.gz
P1__wgts,P1,SB,tumor,dna,fastq,library_id:SB_library;lane:002,/path/to/P1.SB.tumor.dna.wgs.002.R1.fastq.gz;/path/to/P1.SB.tumor.dna.wgs.002.R2.fastq.gz
P1__wgts,P1,SC,tumor,rna,fastq,library_id:SC_library;lane:001,/path/to/P1.SC.tumor.rna.wts.001.R1.fastq.gz;/path/to/P1.SC.tumor.rna.wts.001.R2.fastq.gz
```

### Multiple patients

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
P1__wgts,P1,SA,normal,dna,fastq,library_id:SA_library;lane:001,/path/to/P1.SA.normal.dna.wgs.001.R1.fastq.gz;/path/to/P1.SA.normal.dna.wgs.001.R2.fastq.gz
P1__wgts,P1,SB,tumor,dna,fastq,library_id:SB_library;lane:001,/path/to/P1.SB.tumor.dna.wgs.001.R1.fastq.gz;/path/to/P1.SB.tumor.dna.wgs.001.R2.fastq.gz
P2__wgts,P2,SA,normal,dna,fastq,library_id:SA_library;lane:001,/path/to/P2.SA.normal.dna.wgs.001.R1.fastq.gz;/path/to/P2.SA.normal.dna.wgs.001.R2.fastq.gz
P2__wgts,P2,SB,tumor,dna,fastq,library_id:SB_library;lane:001,/path/to/P2.SB.tumor.dna.wgs.001.R1.fastq.gz;/path/to/P2.SB.tumor.dna.wgs.001.R2.fastq.gz
```

### Column descriptions

| Column        | Description                                                          |
| ------------- | -------------------------------------------------------------------- |
| group_id      | Group ID for a set of samples and inputs                             |
| subject_id    | Subject/patient ID                                                   |
| sample_id     | Sample ID                                                            |
| sample_type   | Sample type: `tumor`, `normal`                                       |
| sequence_type | Sequence type: `dna`, `rna`                                          |
| filetype      | File type: e.g. `fastq`, `bam`, `bai`                                |
| filepath      | Absolute filepath to input file (can be local filepath, URL, S3 URI) |

The identifiers provided in the samplesheet are used to set output file paths:

- `group_id`: top-level output directory for analysis files e.g. `output/COLO829_example/`
- tumor `sample_id`: output prefix for most filenames e.g. `COLO829T.purple.sv.vcf.gz`
- normal `sample_id`: output prefix for some filenames e.g. `COLO829R.cobalt.ratio.pcf`

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision v0.3.1 \
  --mode <wgts|targeted> \
  --genome <GRCh37_hmf|GRCh38_hmf> \
  --input samplesheet.csv \
  --outdir <output_directory>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

> [!NOTE]
> When oncoanalyser is run, it will retrieve all reference data it requires to perform the requested analysis. When
> running oncoanalyser more than once, it is strongly recommended to pre-stage reference data locally to avoid it being
> retrieved multiple times by oncoanalyser. See [Staging reference data](#staging-reference-data).

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/oncoanalyser -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/oncoanalyser
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/oncoanalyser releases page](https://github.com/nf-core/oncoanalyser/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Advanced usage

### Selecting processes

Most of the major components in oncoanalyser can be skipped using `--processes_exclude` (the full list of available
processes can be view [here](https://github.com/nf-core/oncoanalyser/blob/v0.3.1/lib/Constants.groovy#L36-L54)).
Multiple processes can be given as comma-separated list. While there are some use-cases for this feature (e.g. skipping
resource intensive processes such as VIRUSBreakend), it becomes more powerful when combined with existing inputs as
described in the follow section.

> [!WARNING]
> When skipping components, no checks are done to identify orphan processes in the execution DAG or for redundant
> processes.

### Existing inputs

The oncoanalyser pipeline has been designed to allow entry at arbiturary points and is particularly useful in
situtations where previous outputs exist and re-running oncoanalyser is desired (e.g. to subsequently execute an
optional sensor or use an upgrade component such as PURPLE). The primary advantage of this approach is that only the
required processes are executed, which can greatly reduce runtimes by skipping unneccessary processes.

In order to effectively utilise this feature, existing inputs must be set in the [samplesheet](#samplesheet) and the
appropriate [processes selected](#selecting-processes). Take the below example where existing PURPLE inputs are used so
that all upstream variant calling can be skipped:

```csv title='samplesheet.existing_purple.csv'
P1__wgts,P1,SA,normal,dna,bam,/path/to/P1.SA.normal.dna.wgs.bam
P1__wgts,P1,SB,tumor,dna,bam,/path/to/P1.SB.tumor.dna.wgs.bam
P1__wgts,P1,SB,tumor,dna,purple_dir,/path/to/P1.purple_dir/
```

> [!NOTE]
> The original source input file (i.e. BAM or FASTQ) must always be provided for oncoanalyser to infer the correct
> analysis type.

And now run and skip variant calling:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision v0.3.1 \
  --mode wgts \
  --processes_exclude amber,cobalt,gridss,gripss,sage,pave \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

> [!WARNING]
> Providing existing inputs will cause oncoanalyser to skip the corresponding process but _not any_ of the upstream
> processes.

### Configuring reference data

All reference data can be configured as needed. These are defined in various locations:

| Reference data          | Filepath                  | Note                                    |
| ----------------------- | ------------------------- | --------------------------------------- |
| hmftools resource files | `conf/hmf_data.config`    | Paths relative to data bundle directory |
| panel resource files    | `conf/panel_data.config`  | Paths relative to data bundle directory |
| Genomes and indexes     | `conf/hmf_genomes.config` | Absolute paths                          |

To override hmftools resource files (e.g. driver gene panel), [stage the bundle](#staging-reference-data) locally then
copy in the desired file(s) and update `conf/hmf_data.config` accordingly. The local custom bundle must be provided to
oncoanalyser with the `--ref_data_hmf_data_path` CLI option. The same approach is followed for customising panel
resource files, configuring `conf/panel_data.config` and supplying with `--ref_data_panel_data_path` instead.

The path or URI to the VIRUSBreakend database can also be explicitly set with `--ref_data_virusbreakenddb_path`.
Configuring custom genomes uses a different approach to align with the existing concepts in nf-core.

#### Custom genomes

It is strongly recommended to use the Hartwig-distributed reference genomes for alignments
([GRCh37](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/ref_genome/37) or
[GRCh38](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/ref_genome/38)). If there is no
other option than to use a custom genome, one can be configured with the following process:

```text title='genome.custom.config'
params {
    genomes {
        CustomGenome {
            fasta           = "/path/to/CustomGenome/custom_genome.fa"
            fai             = "/path/to/CustomGenome/samtools_index/1.16/custom_genome.fa.fai"
            dict            = "/path/to/CustomGenome/samtools_index/1.16/custom_genome.fa.dict"
            bwa_index       = "/path/to/CustomGenome/bwa_index/0.7.17-r1188/"
            bwa_index_bseq  = "/path/to/CustomGenome/bwa_index/2.2.1/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.0123"
            bwa_index_biidx = "/path/to/CustomGenome/bwa_index/2.2.1/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwt.2bit.64"
            bwa_index_image = "/path/to/CustomGenome/bwa_index_image/0.7.17-r1188/custom_genome.fa.img"
            gridss_index    = "/path/to/CustomGenome/gridss_index/2.13.2/custom_genome.fa.gridsscache"
            star_index      = "/path/to/CustomGenome/star_index/gencode_38/2.7.3a/"
        }
    }
}
```

Run a custom genome with the above configuration and below command

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision v0.3.1 \
  -config genome.custom.config \
  --mode wgts \
  \
  --genome CustomGenome \
  --genome_version <37|38> \
  --genome_type <alt|no_alt> \
  --force_genome \
  \
  --input samplesheet.csv \
  --outdir output/
```

> [!WARNING]
> RNA alignment with STAR must use an index generated from a matching Ensembl release version (GRCh37: v74; GRCh38:
> v104).

#### Staging reference data

Please refer to [REFERENCE_DATA.md](https://github.com/nf-core/oncoanalyser/REFERENCE_DATA.md).

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
