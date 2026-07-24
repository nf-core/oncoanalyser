# FAQ and troubleshooting

- [I only want variant calls, where can I find this data?](#i-only-want-variant-calls-where-can-i-find-this-data)
- [Why does `oncoanalyser` call too many / too few variants than another
  pipeline?](#why-does-oncoanalyser-call-too-many--too-few-variants-than-another-pipeline)
- [My compute environment does not allow Docker](#my-compute-environment-does-not-allow-docker)
- [Running `oncoanalyser` offline](#running-oncoanalyser-offline)
- [Network timeout](#network-timeout)
- [Placing `oncoanalyser` CLI arguments into a configuration
  file](#placing-oncoanalyser-cli-arguments-into-a-configuration-file)
- [Errors and navigating the `work/` directory](#errors-and-navigating-the-work-directory)
- [Resuming runs in Google Batch](#resuming-runs-in-google-batch)

## I only want variant calls, where can I find this data?

Variant calls can be found in the following files:

- PURPLE VCF/TSV files: purity/ploidy adjusted SNV, INDEL, SV, and CNV calls
- SAGE VCF files: raw SNV and INDEL calls
- ESVEE VCF files: raw SV calls

For descriptions of each file, please see the [Output](https://nf-co.re/oncoanalyser/docs/output/) tab.

If you only want to run the variant calling steps, you can either manually select the variant calling processes or
exclude downstream processes (see: [Process selection](./#process-selection)). Using manual process selection for
example, you would run `oncoanalyser` with the below command (assuming starting from FASTQ for DNA sequencing data):

```bash
nextflow run nf-core/oncoanalyser \
  -revision 3.0.0 \
  -profile docker \
  --mode wgts \
  --processes_manual alignment,redux,amber,cobalt,sage,pave,esvee,purple \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

## Why does `oncoanalyser` call too many / too few variants than another pipeline?

The `oncoanalyser` pipeline uses variants with > 2% VAF. Other pipelines may have different assumptions which may cause
differences in samples with low tumor purity or a high number of subclonal variants.

## My compute environment does not allow Docker

Docker is not allowed on some compute environments (especially HPCs) as it runs a daemon as root which is deemed a
security issue. In these cases, using Singularity is recommended by providing `-profile singularity` when running
`oncoanalyser` In these cases, it is recommended to use Singularity (now known as Apptainer) images that are cached for
offline use (see the [Downloading Apptainer
containers](https://nf-co.re/docs/nf-core-tools/pipelines/download#downloading-apptainer-containers)).

:::warning

When manually downloading singularity images, do not execute multiple `singularity pull` commands in parallel. E.g. do
not pull different singularity images in separate terminal sessions on the same compute environment. This will result in
a "[no descriptor found for reference](https://github.com/apptainer/singularity/issues/4555)" error.

:::

:::tip

Docker images can be [pulled with Singularity](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html)
using `singularity pull --name <output_path> <docker_image_url>`.

:::

## Running `oncoanalyser` offline

Sometimes you may need to run `oncoanalyser` on a cloud VM or HPC system with no internet connection. To do this, you
will need to: 1) [manually set up reference data](./#staging-reference-data) and 2) run `oncoanalyser` once (e.g. using
the test profile) to cache Nextflow dependencies and download container images (see nf-core docs [for running pipelines
offline](https://nf-co.re/docs/usage/getting_started/offline)).

Additionally, you may want to add the below item to config your config file:

```groovy
env {
  // Disable automatic update checks, prevents downloading dependencies, execute nextflow using locally available resources
  NXF_OFFLINE = 'true'

  // If NXF_OFFLINE doesn't work, reduce the http timeout to basically zero so that nextflow doesn't hang and throw this error:
  // "Failed to connect to www.nextflow.io port 443 after 300479 ms: Timeout was reached"
  NXF_HTTP_TIMEOUT = '1ms'

  // Nextflow creates and caches dependencies to the '.nextflow/` dir in the current working dir
  NXF_HOME = "/path/to/.nextflow/"
}
```

## Network timeout

The `oncoanalyser` pipeline may time out if pulling containers takes too long. To fix this, increase the network timeout
in the config file (see the [Nextflow config docs](https://www.nextflow.io/docs/latest/reference/config.html) to
configure other container platforms):

```groovy title='timeout.config'
singularity { // If using singularity
  pullTimeout = '2h'
}

docker { // If using Docker
  pullTimeout = '2h'
}

env { // Shell environment variables
  NXF_HTTP_TIMEOUT = '2h'
}
```

Network timeout may also occur when downloading reference data. While the above solution might also work, we recommend
downloading and setting up reference data [manually](./#staging-reference-data) instead.

## Placing `oncoanalyser` CLI arguments into a configuration file

Almost all `oncoanalyser` arguments in the [Parameters tab](https://nf-co.re/oncoanalyser/parameters/) can be placed
in a config file.

For example, the `oncoanalyser` arguments which start with `--` in this command:

```shell
nextflow run nf-core/oncoanalyser \
  -revision 3.0.0 \
  -profile docker \
  -config refdata.config \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

can be specified in a config file by stripping the `--` like so:

```groovy title='params.config'
params {
    mode   = "wgts"
    genome = "GRCh38_hmf"
    input  = "samplesheet.csv"
    outdir = "outdir/"
}
```

and provided as a config file when running `oncoanalyser`:

```shell
nextflow run nf-core/oncoanalyser \
  -revision 3.0.0 \
  -config refdata.config \
  -config params.config \
  -profile docker \
  <...>
```

:::tip

The `-config` Nextflow argument can be used multiple times to provide multiple config file.

:::

## Errors and navigating the `work/` directory

When `oncoanalyser` crashes, you may need to further investigate error messages in the `.nextflow.log` files or the
`work` directory.

The `work` directory contains the run scripts, logs, and input/output files for each process. Once the process is done
running, only the output files are 'published' (copied) to the final output directory (as specified by `--outdir`).

Below is an example `work` directory for one process. Error messages will typically be found in `.command.log`,
`.command.err` or `.command.out` log files. You can send these logs / error messages for example to the `oncoanalyser`
[Slack channel](https://nfcore.slack.com/channels/oncoanalyser), as an issue on the
[oncoanalyser](https://github.com/nf-core/oncoanalyser) or [WiGiTS](https://github.com/hartwigmedical/hmftools) GitHub
repositories.

```shell
work/
├── e5
│   └── f6e2e8f18ef70add9349164d5fb37e
│       ├── .command.sh     # Bash script used to run the process *within the container*
│       ├── .command.run    # Bash script used to run the process in the host machine
│       ├── .command.begin
│       ├── .command.log    # All log messages (combination of stdout and stderr). Might not exist for some executors
│       ├── .command.err    # stderr log messages
│       ├── .command.out    # stdout log messages
│       ├── .command.trace  # Compute resource usage stats
│       ├── .exitcode       # Exit code
│       ├── <...>           # Input/output files or directories
│       └── versions.yml    # WiGiTS tool version
<...>
```

The `work/` directory can be hard to navigate due to the `<short_hash>/<long_hash>` structure. These hashes are shown
(truncated) in the console while running `oncoanalyser` (but can also be found in the `.nextflow.log` files):

```shell
[e5/f6e2e8] process > NFCORE_ONCOANALYSER:WGTS:REDUX_PROCESSING:REDUX (<group_id>_<sample_id>)     [100%] 2 of 2 ✔
```

Otherwise, you can use a utility like [tree](<https://en.wikipedia.org/wiki/Tree_(command)>) to show the directory
structure, which allows you to manually find the target process directory.

```shell
outdir/coloMini/logs/NFCORE_ONCOANALYSER:WGTS:REDUX_PROCESSING:REDUX.coloMini_coloMiniT.command.log
```

## Resuming runs in Google Batch?

When resuming with runs in Google Batch (using `-resume`), you will need to enable overwriting of the `pipeline_info`
files (performance and run time stats) as shown below. By default, these files are not overwritten thus preventing
`oncoanalyser` from starting.

```groovy
timeline {
    enabled   = true
    overwrite = true
    file      = "${params.outdir}/pipeline_info/execution_timeline.html"
}

report {
    enabled   = true
    overwrite = true
    file      = "${params.outdir}/pipeline_info/execution_report.html"
}

trace {
    enabled   = true
    overwrite = true
    file      = "${params.outdir}/pipeline_info/execution_trace.txt"
}

dag {
    enabled   = true
    overwrite = true
    file      = "${params.outdir}/pipeline_info/pipeline_dag.svg"
}
```
