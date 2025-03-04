# nf-core/oncoanalyser: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/oncoanalyser/usage](https://nf-co.re/oncoanalyser/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Table of contents

**TEMPORARY: REMOVE THIS BEFORE PUBLISHING TO NF-CORE**

<!-- TOC -->
* [nf-core/oncoanalyser: Usage](#nf-coreoncoanalyser--usage)
  * [:warning: Please read this documentation on the nf-core website: https://nf-co.re/oncoanalyser/usage](#-warning--please-read-this-documentation-on-the-nf-core-website--httpsnf-coreoncoanalyserusage)
  * [Table of contents](#table-of-contents)
  * [Introduction](#introduction)
  * [Getting started](#getting-started)
  * [Running the pipeline](#running-the-pipeline)
    * [Command line interface (CLI)](#command-line-interface--cli-)
    * [Outputs](#outputs)
    * [Reusing CLI arguments](#reusing-cli-arguments)
    * [Versions and reproducibility](#versions-and-reproducibility)
  * [Samplesheet](#samplesheet)
    * [Input starting points](#input-starting-points)
      * [FASTQ](#fastq)
      * [BAM / CRAM](#bam--cram)
      * [REDUX BAM](#redux-bam)
    * [Sample setups](#sample-setups)
      * [Tumor/normal DNA](#tumornormal-dna)
      * [Tumor-only DNA](#tumor-only-dna)
      * [Tumor-only RNA](#tumor-only-rna)
      * [Tumor/normal DNA and tumor-only RNA](#tumornormal-dna-and-tumor-only-rna)
      * [Tumor/normal DNA with donor sample](#tumornormal-dna-with-donor-sample)
    * [Multiple samples](#multiple-samples)
  * [Reference data](#reference-data)
    * [Automatic staging of reference data](#automatic-staging-of-reference-data)
    * [Manually staging reference data](#manually-staging-reference-data)
    * [Reference data URLs](#reference-data-urls)
    * [Custom panel reference data](#custom-panel-reference-data)
    * [Custom genomes](#custom-genomes)
  * [Process selection](#process-selection)
    * [Excluding processes](#excluding-processes)
    * [Manual process selection](#manual-process-selection)
    * [Starting from existing inputs](#starting-from-existing-inputs)
  * [Core Nextflow arguments](#core-nextflow-arguments)
    * [`-profile`](#-profile)
    * [`-resume`](#-resume)
    * [`-c`](#-c)
  * [Custom configuration](#custom-configuration)
    * [Compute resources](#compute-resources)
    * [Container images](#container-images)
      * [Setting up Singularity](#setting-up-singularity)
      * [Container image sources](#container-image-sources)
    * [Custom tool arguments](#custom-tool-arguments)
    * [UMI processing](#umi-processing)
    * [nf-core/configs](#nf-coreconfigs)
  * [FAQ and troubleshooting](#faq-and-troubleshooting)
    * [How to start from CRAM?](#how-to-start-from-cram)
    * [How to handle UMIs?](#how-to-handle-umis)
    * [How to use oncoanalyser with a custom panel or whole exome?](#how-to-use-oncoanalyser-with-a-custom-panel-or-whole-exome)
    * [Incompatible BAM files](#incompatible-bam-files)
    * [I want to store the output BAMs. Why are there only REDUX BAM(s) with additional files?](#i-want-to-store-the-output-bams-why-are-there-only-redux-bam--s--with-additional-files)
    * [I only want variant calls](#i-only-want-variant-calls)
    * [Why does `oncoanalyser` call too many / too few variants than another pipeline?](#why-does-oncoanalyser-call-too-many--too-few-variants-than-another-pipeline)
    * [My HPC does not allow Docker](#my-hpc-does-not-allow-docker)
    * [Firewall prevents `oncoanalyser` from pulling Singularity images](#firewall-prevents-oncoanalyser-from-pulling-singularity-images)
    * [Network timeout](#network-timeout)
    * [Run fails due to insufficient CPUs/RAM/disk](#run-fails-due-to-insufficient-cpusramdisk)
    * [Can `oncoanalyser` CLI arguments be put in a config file?](#can-oncoanalyser-cli-arguments-be-put-in-a-config-file)
    * [Errors and navigating the `work/` directory](#errors-and-navigating-the-work-directory)
  * [Azure resource requests](#azure-resource-requests)
  * [Running in the background](#running-in-the-background)
  * [Nextflow memory requirements](#nextflow-memory-requirements)
<!-- TOC -->

## Introduction

The `oncoanalyser` pipeline typically runs from  **FASTQ**, **BAM** or **CRAM** [input files](#input-starting-points), supports most
**GRCh37** and **GRCh38** human reference genome builds, and supports **UMI** ([unique molecular identifier](#umi-processing)) processing
for DNA sequencing data.

The pipeline supports **WGTS** (whole genome and/or transcriptome) and **targeted** panel workflow modes; and supports **tumor/normal**
(with optional donor sample) and **tumor-only** [sample setups](#sample-setups). The below table summarises the supported analyses:

| Sample setup           | WGTS workflow  | Targeted workflow |
|------------------------|----------------|-------------------|
| Tumor / Normal         | DNA and/or RNA | DNA only          |
| Tumor / Normal / Donor | DNA and/or RNA | DNA only          |
| Tumor-only             | DNA and/or RNA | DNA and/or RNA    |

## Getting started

Assuming that [**Nextflow**](https://www.nextflow.io/docs/latest/install.html) and
[**Docker**](https://docs.docker.com/engine/install/) or [**Singularity**](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#)
are installed, setting up and running `oncoanalyser` involves the following steps:

1. Download and configure [**reference data**](#reference-data) (e.g. reference genome)
2. (Optional) [**Other configuration**](#custom-configuration) (e.g. compute resources)
3. Create [**samplesheet**](#samplesheet) specifying input files
4. [**Run**](#running-the-pipeline) `oncoanalyser`

:::tip

Jump to [**FAQ and troubleshooting**](#faq-and-troubleshooting).

:::

## Running the pipeline

### Command line interface (CLI)
A typical command for running `oncoanalyser` is shown below:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/ \
  -config oncoanalyser.config ## Optional; see below note
```

`oncoanalyser` supports two analysis `--mode`s:
- `wgts`: whole genome and/or transcriptome sequencing data
- `targeted`: targeted sequencing data (panels or exomes)

:::tip

**Nextflow** arguments have **one hyphen** (`-`) and are detailed in the [**Nextflow CLI documentation**](https://www.nextflow.io/docs/latest/reference/cli.html).

**Oncoanalyser** arguments have **two hyphens** (`--`) and are detailed in the [**Parameters tab**](https://nf-co.re/oncoanalyser/parameters).

:::

:::note

It is strongly recommended to locally stage and configure reference data paths (e.g. reference genome) with `-config <file>`
(see: [**Manually setting up reference data**](#manually-staging-reference-data)). If not configured, `oncoanalyser` will automatically
re-download reference data for every run executed from a different working directory, which can lead to unnecessary disk/network usage.

:::

### Outputs

`oncoanalyser` will create the following files in your working directory:

```bash
work           # Directory containing the nextflow working files
<OUTDIR>       # Finished results in specified location (defined with --outdir)
.nextflow_log  # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

Descriptions of each file in `output/` is described in the [**Output**](https://nf-co.re/oncoanalyser/docs/output/) documentation.

### Reusing CLI arguments

To use the same CLI arguments across multiple runs, it may be handy to specify these in a `yaml` or `json` file via `-params-file <file>`.
The [above command](#command-line-interface--cli-) would have the equivalent `yaml` file:

```yaml title="params.yaml
mode:   'wgts'
genome: 'GRCh38_hmf'
input:  './samplesheet.csv'
outdir: './output/'
<...>
```

and be run using this command:

```bash
nextflow run nf-core/oncoanalyser -profile docker -params-file params.yaml
```

You can also generate such `yaml`/`json` files via [nf-core/launch](https://nf-co.re/launch).

### Versions and reproducibility

It is recommended to specify a version/tag with e.g. `-revision 2.0.0` when running `oncoanalyser`. This ensures the same pipeline version
will run even if there have been changes to the code since. Tags can be found on the [nf-core/oncoanalyser releases page](https://github.com/nf-core/oncoanalyser/tags).
The pipeline version will be logged in `<outdir>/pipeline_info/software_versions.yml`. [Reusing parameter files](#reusing-cli-arguments) can
also further assist in reproducibility.

When running `nextflow run nf-core/oncoanalyser <args>` without `-r`/`-revision`, Nextflow automatically pulls and locally caches the
latest `oncoanalyser` version from [GitHub](https://github.com/nf-core/oncoanalyser). Subsequent runs will use the cached version if
available, even if the pipeline has been updated since. To ensure you're running the latest `oncoanalyser` version, you can manually update
the cached code like so:

```bash
nextflow pull nf-core/oncoanalyser
```

## Samplesheet

The samplesheet is a comma separated table with the following columns:

| Column          | Description                                                                                                                                                                                                                                  |
|:----------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `group_id`      | Groups `sample_id` entries (e.g. tumor DNA, normal DNA, tumor RNA for one patient) into the same analysis.</br>Output files are also organized by `group_id` (see [Output](https://nf-co.re/oncoanalyser/1.0.0/docs/output/) documentation). |
| `subject_id`    | Must correlate exactly with `group_id`                                                                                                                                                                                                       |
| `sample_id`     | Sample ID                                                                                                                                                                                                                                    |
| `sample_type`   | Can be: `tumor` or `normal`                                                                                                                                                                                                                  |
| `sequence_type` | Can be: `dna` or `rna`                                                                                                                                                                                                                       |
| `filetype`      | Use `bam`/`bai` for BAMs or CRAMs.<br/>Use `fastq` for FASTQs.</br>The full list of valid values can be viewed [here](http://github.com/nf-core/oncoanalyser/blob/master/lib/Constants.groovy#L58)                                           |
| `filepath`      | Absolute filepath to input file. Can be local file path, URL, or S3 URI                                                                                                                                                                      |
| `info`          | Sequencing library and lane info for [FASTQ](#fastq) files. Column does not need to be provided if not running from FASTQ                                                                                                                    |

:::note

`subject_id` is reserved for implementing more complex sample setups in the future, but currently has no functional purpose.
For now, `subject_id` must correlate exactly with `group_id`. E.g. if `group_id` is `PATIENT1_WGS` and `subject_id` is `PATIENT1`,
this must be the case for every row in the samplesheet.

:::

### Input starting points

Below subsections provide example samplesheets for starting from FASTQ, BAM, CRAM, or REDUX BAM.

#### FASTQ

To run from FASTQ, provide:
- In column `info`: the sequencing library and lane info separated by `;`
- In column `filepath`, the forward ('R1') and reverse ('R2') FASTQ files separated by `;`

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath,info
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,fastq,/path/to/PATIENT1-T_S1_L001_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L001_R2_001.fastq.gz,library_id:S1;lane:001
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,fastq,/path/to/PATIENT1-T_S1_L002_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L002_R2_001.fastq.gz,library_id:S1;lane:002
```

:::note

Currently only gzip compressed, non-interleaved pair-end FASTQ files are currently supported

:::

#### BAM / CRAM

To run from BAM, specify `bam` in column `filetype`:

```csv title="samplesheet.bam.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
```

BAM indexes (.bai files) are expected to be in the same directory as the BAM files. You can also explicitly provide the BAM index path by
specifying `bai` in column `filetype`:

```csv title="samplesheet.bam_bai.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/other/dir/PATIENT1-T.dna.bam.bai
```

To run from CRAM, simply provide the CRAM and/or CRAM index to `bam`/`bai` ibn column `filetype`:
```csv title="samplesheet.cram.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.cram
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/other/dir/PATIENT1-T.dna.cram.crai
```

#### REDUX BAM
For DNA sequencing analyses, read alignment with [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) and read pre-processing with
[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) are the pipeline steps that take the most time and compute resources.

We can re-run `oncoanalyser` from a REDUX BAM if it already exists by specifying `bam_redux` in column `filetype`. We can optionally
explicitly provide the BAM index with `filetype` as `bai`.

```csv title="samplesheet.redux_bam.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam_redux,/path/to/PATIENT1-T.dna.redux.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/other/dir/PATIENT1-T.dna.redux.bam.bai
```

The `*.jitter_params.tsv` and `*.ms_table.tsv.gz` REDUX output files are expected to be in the same directory as the REDUX BAM, and are
required by [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage). If these files are located elsewhere, their paths can
be explicitly provided by specifying `redux_jitter_tsv` and `redux_ms_tsv` under `filetype`:

```csv title="samplesheet.redux_inputs.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam_redux,/path/to/PATIENT1-T.dna.redux.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,redux_jitter_tsv,/other/dir/PATIENT1-T.dna.jitter_params.tsv
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,redux_ms_tsv,/path/dir/PATIENT1-T.dna.ms_table.tsv.gz
```

### Sample setups

Providing `sample_type` and `sequence_type` in different combinations allows `oncoanalyser` to run in different sample setups. The below
samplesheet examples use BAM files, but different sample setups can also be specified for FASTQ or CRAM files.

#### Tumor/normal DNA
```csv title="samplesheet.tn_dna.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
```

#### Tumor-only DNA

```csv title="samplesheet.to_dna.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
```

#### Tumor-only RNA
```csv title="samplesheet.to_rna.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.bam
```

#### Tumor/normal DNA and tumor-only RNA

```csv title="samplesheet.wgts.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.bam
```

#### Tumor/normal DNA with donor sample

Adding a donor sample is useful when a patient has had e.g. a bone marrow transplant. This allows for germline variant subtraction from both
the patient's normal sample and the bone marrow donor's normal sample (typical blood in both cases).

To add a donor sample, specify `donor` in column `sample_type` and a different sample ID in column `sample_id`:

```csv title="samplesheet.tn_with_donor.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
PATIENT1,PATIENT1,PATIENT1-D,donor,dna,bam,/path/to/PATIENT1-D.dna.bam
```

### Multiple samples

To run with multiple samples, specify a different `group_id`/`subject_id` for each sample:

```csv title="samplesheet.batch.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,bam,/path/to/PATIENT2-T.dna.bam
PATIENT2,PATIENT2,PATIENT2-R,normal,dna,bam,/path/to/PATIENT2-R.dna.bam
```

:::warning

Running a large batch of samples can lead to stalling due to sample runs competing for compute resources. It is therefore recommended to
run each sample in the batch in e.g. a separate cloud VM.

:::

## Reference data

Reference data refers to reference genomes + indexes, [WiGiTS](https://github.com/hartwigmedical/hmftools) resources files, and panel
specific resource files. Descriptions for each file can be found on the
[WiGiTS resource file documentation](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_RESOURCES.md).

If reference data is not configured using `-config <file>`, `oncoanalyser` will automatically download/extract reference data for every
run from a different working directory. Therefore, it is strongly recommended to [manually stage reference data](#manually-staging-reference-data),
especially if you are planning to do multiple `oncoanalyser` runs and/or with different experimental setups (e.g. WGS vs panel).

### Automatic staging of reference data
`oncoanalyser` upon first run will automatically stage (download, extract and cache) default reference data based on the provided samplesheet
and CLI arguments.

For example the below samplehseet and command will stage the `GRCh38_hmf` genome, DNA alignment indexes, and
[WiGiTS](https://github.com/hartwigmedical/hmftools) resources files (but not RNA aligment indexes or panel specific data). Then, the full
pipeline will run.

```csv title="samplesheet.tn_dna.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.bam
```

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

To request `oncoanalyser` to stage reference data without running the full pipeline, add the `--prepare_reference_only` argument:
```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/ \
  --prepare_reference_only
```

Once the above command completes, the prepared reference files can be found in `output/reference_data/2.0.0/<datetimestamp>/`.
The Nextflow `work/` directory can be removed after staging data to free disk space.

### Manually staging reference data

To manually stage reference data, download the files from the link in [**Reference data URLs**](#reference-data-urls).
TAR files will need to be extracted, e.g. using `tar -xzvf <file>`.

Then, create a config file can be used to point to the file paths:

```groovy title="refdata.config"
params {

    genomes {
        GRCh38_hmf {
            fasta            = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
            fai              = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
            dict             = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict"
            img              = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img"

            // Only required when aligning DNA reads from FASTQ
            bwamem2_index    = "/path/to/bwa-mem2_index_dir/"

            // Only required when running Virusbreakend/Virusinterpreter
            gridss_index     = "/path/to/gridss_index_dir/"

            // Only required when aligning RNA reads from FASTQ
            star_index       = "/path/to/star_index_dir/"
        }
    }

    ref_data_hmf_data_path   = "/path/to/hmftools_data_dir/"
}
```

Lastly, provide the config file using `-config <file>` when running `oncoanalyser`:

```bash
nextflow run nf-core/oncoanalyser \
  <...>
  -config refdata.config
```

### Reference data URLs

_GRCh37 genome (Hartwig): `GRCh37_hmf`_

| Type                 | Link                                                                                                                                                                                                |
|:---------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FASTA                | [Homo_sapiens.GRCh37.GATK.illumina.fasta](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/Homo_sapiens.GRCh37.GATK.illumina.fasta)                                      |
| FASTA index          | [Homo_sapiens.GRCh37.GATK.illumina.fasta.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai)          |
| FASTA dictionary     | [Homo_sapiens.GRCh37.GATK.illumina.fasta.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict)        |
| BWA-MEM2 index image | [Homo_sapiens.GRCh37.GATK.illumina.fasta.img](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/bwa_index_image/0.7.17-r1188/Homo_sapiens.GRCh37.GATK.illumina.fasta.img) |
| BWA-MEM2 indexes     | [bwa-mem2_index/2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                              |
| GRIDSS indexes       | [gridss_index/2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                |
| STAR indexes         | [star_index/gencode_19/2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/star_index/gencode_19/2.7.3a.tar.gz)                                              |
| WiGiTS data          | [hmf_pipeline_resources.37_v2.0.0--2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.0.0--2.tar.gz)                            |
| TSO500 panel data    | [hmf_panel_resources.tso500.37_v2.0.0--2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.37_v2.0.0--2.tar.gz)                      |

_GRCh38 genome (Hartwig): `GRCh38_hmf`_

| Type                 | Link                                                                                                                                                                                                                |
|:---------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FASTA                | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna)                                      |
| FASTA index          | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai)          |
| FASTA dictionary     | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict)        |
| BWA-MEM2 index image | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/bwa_index_image/0.7.17-r1188/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.img) |
| BWA-MEM2 indexes     | [bwa-mem2_index/2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                                              |
| GRIDSS indexes       | [gridss_index/2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                                |
| STAR indexes         | [star_index/gencode_38/2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/star_index/gencode_38/2.7.3a.tar.gz)                                                              |
| WiGiTS data          | [hmf_pipeline_resources.38_v2.0.0--2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.0.0--2.tar.gz)                                            |
| TSO500 panel data    | [hmf_panel_resources.tso500.38_v2.0.0--2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.38_v2.0.0--2.tar.gz)                                      |

### Custom panel reference data

To use custom panels, reference data must first be generated using a training procedure detailed
[here](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md). This procedure allows for normalisation of
copy number, tumor mutational burden, and RNA transcript counts, as well as to filter out panel specific artefacts.

:::note

The panel reference data training will soon be integrated into `oncoanalyser` which will greatly simplify the process.

:::

The panel reference data paths must then be provided in a config file:

```groovy title="panel.config"
params {
    //
    ref_data_panel_data_path = "/path/to/panel_resources/"

    // These are relative paths within the dir provided by `ref_data_panel_data_path`
    panel_data_paths {

        custom_panel { // This is the name that should be passed to the `--panel` argument

            // Can be '37' or '38'
            '38' {

                driver_gene_panel           = 'common/DriverGenePanel.custom_panel.38.tsv'
                sage_actionable_panel       = 'variants/ActionableCodingPanel.custom_panel.38.bed.gz'
                sage_coverage_panel         = 'variants/CoverageCodingPanel.custom_panel.38.bed.gz'
                pon_artefacts               = 'variants/pon_artefacts.custom_panel.38.tsv.gz'
                target_region_bed           = 'copy_number/target_regions_definition.custom_panel.38.bed.gz'
                target_region_normalisation = 'copy_number/cobalt_normalisation.custom_panel.38.tsv'
                target_region_ratios        = 'copy_number/target_regions_ratios.custom_panel.38.tsv'
                target_region_msi_indels    = 'copy_number/target_regions_msi_indels.custom_panel.38.tsv'

                // Optional. If no RNA in panel, these can be omitted by providing in empty list, e.g.:
                // isofox_tpm_norm = []
                isofox_tpm_norm             = 'rna_resources/isofox.gene_normalisation.custom_panel.38.csv'
                isofox_gene_ids             = 'rna_resources/custom_panel.rna_gene_ids.csv'
                isofox_counts               = 'rna_resources/read_93_exp_counts.38.csv'
                isofox_gc_ratios            = 'rna_resources/read_93_exp_gc_ratios.38.csv'
            }
        }
    }
}
```

Lastly, run `oncoanalyser`:
- Provide both the general and panel reference data config files to `-config`
- Pass the panel name to `--panel`. This should match the name defined in the panel resources config file
- Provide argument `--force_panel` if `--panel` is not `tso500` (this is currently the only supported panel type)

```bash
nextflow run nf-core/oncoanalyser \
  --panel custom_panel \
  --force_panel \
  -config refdata.config \
  -config panel.config \
  --mode targeted \
  <...>
```

### Custom genomes

It is strongly recommended to use a Hartwig-distributed reference genome for alignments and subsequent analysis
(`GRCh37_hmf` or `GRCh38_hmf`). Where it is not feasible to do so, a custom genome can instead be used by providing the
relevant FASTA file in a configuration file:

```groovy title='genome.custom.config'
params {
    genomes {
        CustomGenome {
            fasta = "/path/to/custom_genome.fa"
        }
    }
}
```

Each index required for the analysis will first be created before running the rest of `oncoanalyser` with the following command:

:::warning

For GRCh38 genome builds, It is strongly recommended to remove or mask HLA class I alt contigs. These lead to poor variant calling in HLA
class I genes and poor HLA typing. We will soon release a standalone utility to remap HLA contigs back to the original BAM for users who
want to start from BAM.

:::

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
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

:::tip

You can first generate the required indexes by setting `--prepare_reference_only` and then provide the prepared reference files to
`oncoanalyser` through a custom config file. This avoids having to regenerate indexes for each new analysis.

:::

Creation of a STAR index also requires transcript annotations, please provide either of the following GTF files via the
`--ref_data_genome_gtf` option after decompressing:

- GRCh37: [GENCODE v37 (Ensembl v74) annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz)
- GRCh38: [GENCODE v38 (Ensembl v104) annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)

:::warning

STAR index must use transcript annotations from Ensembl versions that match WiGiTS resource data (GRCh37: v74; GRCh38: v104).

:::

When creating indexes for reference genomes with alternative haplotypes, an ALT file must be given with
`--ref_data_genome_alt`. Importantly, a STAR index will not be generated for reference genomes with alternative
haplotypes since this requires careful processing and is hence left to the user.

## Process selection

It is possible to exclude or include specific processes when running `oncoanalyser`.  The full list of processes that can be selected
viewed [here](https://github.com/nf-core/oncoanalyser/blob/master/lib/Constants.groovy#L36)).

### Excluding processes

Most of the major components in `oncoanalyser` can be skipped using `--processes_exclude`. Multiple processes can be given as a
comma-separated list.

You may want to do this to skip resource intensive processes such as VIRUSBreakend, or you simply have no use for e.g. the ORANGE report.
The command to run `oncoanalyser` would then be:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/ \
  --processes_exclude virusinterpreter,orange
```

:::warning

When skipping components no checks are done to identify orphan processes in the execution DAG or for redundant processes.

:::

### Manual process selection

We can use argument `--processes_manual` to enable manual process selection, and use argument `--processes_include <process_1,process_2>` to
choose which processes to run.

One use case would be to run processes which are turned off by default, such as neoepitope calling with
[NEO](https://github.com/hartwigmedical/hmftools/tree/master/neo). To do this, provide the below example samplesheet:

```csv title='samplesheet.neo_inputs.csv'
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.wgs.bam
```

Then, run `oncoanalyser` with the `neo` process selected as well as all required upstream processes:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.neo_inputs.csv \
  --outdir output/ \
  --processes_manual \
  --processes_include neo,redux,amber,cobalt,sage,pave,esvee,purple,lilac,isofox
```

:::warning

It is the user's responsibility to select the required upstream processes for a downstream process to run. If not all required processes
are not selected, `oncoanalyser` will not crash but will simply finish without the downstream process running.

:::

### Starting from existing inputs
It is possible to run downstream tools given the outputs from upstream `oncoanalyser` processes. For example, you want to run
[NEO](https://github.com/hartwigmedical/hmftools/tree/master/neo) given
existing outputs from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple),
[LILAC](https://github.com/hartwigmedical/hmftools/tree/master/lilac) and
[ISOFOX](https://github.com/hartwigmedical/hmftools/tree/master/isofox). To do this, provide the below example samplesheet:

```csv title='samplesheet.neo_inputs.csv'
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-R,normal,dna,bam,/path/to/PATIENT1-R.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-RNA,tumor,rna,isofox_dir,/path/to/PATIENT1.isofox_dir/
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,purple_dir,/path/to/PATIENT1.purple_dir/
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,lilac_dir,/path/to/PATIENT1.lilac_dir/
```

:::note

The original source input file (i.e. BAM or FASTQ) must be provided for `oncoanalyser` to infer the correct analysis type.

:::

Then, run `oncoanalyser`, skipping variant calling by using `--processes_exclude` (see [Excluding process](#excluding-processes)):

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.existing_purple.csv \
  --outdir output/ \
  --processes_exclude redux,amber,cobalt,esvee,sage,pave
```

:::warning

Providing existing inputs will cause `oncoanalyser` to skip the corresponding process but _not any_ of the upstream
processes. It is the responsibility of the user to skip all relevant processes.

:::

## Core Nextflow arguments

:::note

These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker,
Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info

We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible,
Conda is also supported.

:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it
runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your
system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is
_not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

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

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing
from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well.
For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the
[nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

Custom configuration can be provided to `oncoanalyser` by providing a config file to the CLI argument `-config <file>` or `-c <file>`.
Syntax and examples of config items are described in the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and
[nf-core documentation](https://nf-co.re/usage/configuration). Below a subsections describe common use cases for custom configuration.

:::tip

See [this example config file](https://github.com/nf-core/oncoanalyser/blob/update-docs-v2.0/conf/hartwig.config) for common configuration
options.

:::

### Recommended configuration

Recommended configurations for oncoanalyser can readily be loaded using the ```-profile hartwig``` command line argument, whenever launching a run.
Those specifications can currently be found [here](https://github.com/nf-core/oncoanalyser/blob/update-docs-v2.0/conf/hartwig.config) (will be replaced with nf-core repo path)


### Compute resources

Depending on the compute environment/setup, the default compute resources (e.g. CPUs, RAM, disk space) configured in `oncoanalyser` may not
be sufficient for one or more processes. To change the resource requests, please see the
[max resources](https://nf-co.re/docs/usage/configuration#max-resources) and
[tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) sections of the nf-core website.

Note that for most of the steps in the pipeline, failing jobs with any of the error codes specified in
[this](https://github.com/nf-core/oncoanalyser/blob/master/conf/base.config) base config file will automatically be resubmitted with higher
requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

### Container images

You may want to change which container or conda environment uses for a particular process (e.g. due to a newer tool
version being available). Please see [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) for
instructions.

#### Setting up Singularity

Some restricted compute environments (e.g. HPC) do not allow Docker to be used and/or dynamic download of images upon running
`oncoanalyser`. In these cases, it is recommended to use Singularity (now known as Apptainer) images that are cached for offline use. See
the [Downloading Apptainer containers](https://nf-co.re/docs/nf-core-tools/pipelines/download#downloading-apptainer-containers) section
of the nf-core-tools documentation for details.

:::warning

When manually downloading singularity images, do not execute multiple `singularity pull` commands in parallel. E.g. do not pull different
singularity images in separate terminal sessions on the same compute environment. This will result in a
"[no descriptor found for reference](https://github.com/apptainer/singularity/issues/4555)" error.

:::

:::tip

Docker images can be [pulled with Singularity](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html) using
`singularity pull --name <output_path> <docker_image_url>`

:::

#### Container image sources

Oncoanalyser by default runs each tool using [Docker](https://docs.docker.com/engine/install/) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#) container images which are built by the
[Bioconda recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes) CI/CD infrastructure.
Below is information on how these images can be accessed should you want to download images manually
(e.g. to [set up Singularity for offline use](#setting-up-singularity)).

**Docker (Bioconda)**
- Host: [quay.io](https://quay.io/organization/biocontainers)
- Repo URL example: https://quay.io/repository/biocontainers/hmftools-redux?tab=tags
- Image URI example: `quay.io/biocontainers/hmftools-redux:1.1--hdfd78af_1`

**Singularity (Bioconda)**
- Host: [Galaxy Project](https://depot.galaxyproject.org/singularity/)
- Image URL example: https://depot.galaxyproject.org/singularity/hmftools-redux:1.1--hdfd78af_1

**Bioconda recipes** for the above containers are found here:
- Host: https://github.com/bioconda/bioconda-recipes/
- Recipe example: https://github.com/bioconda/bioconda-recipes/tree/master/recipes/hmftools-redux

**Docker** images built by Hartwig's CI/CD infrastructure are also available, intended for beta releases and not used by default in `oncoanalyser`
- Host: [Dockerhub](https://hub.docker.com/r/hartwigmedicalfoundation)
- Repo URL example: https://hub.docker.com/r/hartwigmedicalfoundation/redux/tags
- Image URI example: `docker.io/hartwigmedicalfoundation/sage:1.1`

### Custom tool arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines
provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the
[customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### UMI processing

Unique molecular identifiers (UMI) allow for read deduplication and error correction. UMI processing is performed by
[fastp](https://github.com/OpenGene/fastp?tab=readme-ov-file#unique-molecular-identifier-umi-processing) for FASTQ files, and
[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux#deduplication) for BAM files. Depending on the presence/format of your
UMI strings, you may need to configure one or more of these arguments:

```groovy title='umi.config'
params {
    // For FASTQ files
    fastp_umi = true                // Enable UMI processing by fastp
    fastp_umi_location = "per_read" // --umi_loc fastp arg
    fastp_umi_length = 7            // --umi_len fastp arg
    fastp_umi_skip = 0              // --umi_skip fastp arg

    // For BAM files
    redux_umi = true                // Enable UMI processing by REDUX
    redux_umi_duplex_delim = "_"    // Duplex UMI delimiter
}
```

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be
running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config
file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your
pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of
your config file, associated documentation file (see examples in
[`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending
[`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own
configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the
[`#configs` channel](https://nfcore.slack.com/channels/configs).


## FAQ and troubleshooting

### How to start from CRAM?
Simply specify a CRAM path instead of a BAM path in the sample sheet. See section [Input starting points: BAM / CRAM](#bam--cram).

### How to handle UMIs?
UMI processing can be enabled and configured via a config file. See section [UMI processing](#umi-processing).

### How to use oncoanalyser with a custom panel or whole exome?
`oncoanalyser` currently has built-in support for the TSO500 panel. For custom panels however, reference data must first be generated using
a training procedure detailed [here](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md). This procedure
allows for normalisation of copy number, TMB, and TPM data as well as filtering of panel specific artefacts.

:::note

We realise that the panel training procedure is not very straight forward and will therefore be integrated into `oncoanalyser` in the next minor release!

:::

### Incompatible BAM files
`oncoanalyser` has been validated on with BAMs aligned with BWA-MEM, BWA-MEM2 and DRAGEN. BAM files from other aligners / sources may be
incompatible with `oncoanalyser` can cause the pipeline to crash.

One requirement for example that the mate cigar attribute must be present for any BAM records with paired reads. Non-compatible BAMs may be
rectified using tools such as [Picard FixMateInformation](https://gatk.broadinstitute.org/hc/en-us/articles/360036713471-FixMateInformation-Picard)
routine.

In other cases, converting from BAM back to FASTQ may be required to run `oncoanalyser`.

### I want to store the output BAMs. Why are there only REDUX BAM(s) with additional files?

[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) performs some important read post-processing steps:
- Unmapping of reads in pre-defined problematic regions (extremely high depth, reads often discordant or have long soft clipping). This is
done to remove obvious poor alignments from the BAM prior to running downstream tools
reads are retained in the BAM
- [Read deduplication](https://github.com/hartwigmedical/hmftools/tree/master/redux#deduplication) to form of consensus read with consensus
sequence / base qualities
- Measure the rate of microsatellite errors (see: [jitter modeling](https://github.com/hartwigmedical/hmftools/tree/master/redux#microsatellite-jitter-modelling))
which are store in lookup files (`*.jitter_params.tsv` and `*.ms_table.tsv.gz`)to be used downstream by
[SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage#key-concepts-in-sage) for error-calibrated small variant calling.

It was therefore a choice to provide the user the REDUX BAM (plus TSV files) as output, rather than BAMs from BWA-MEM2 which have
potentially more poor alignments and read duplicates.

:::note

When storing REDUX BAMs, the `*.jitter_params.tsv` and `*.ms_table.tsv.gz` must also be stored!

:::

### I only want variant calls

Variant calls can be found in the following files:
- PURPLE VCF/TSV files: purity/ploidy adjusted SNV, INDEL, SV, and CNV calls
- SAGE VCF files: raw SNV and INDEL calls
- ESVEE VCF files: raw SV calls

For descriptions of each file, please see the [Output](https://nf-co.re/oncoanalyser/docs/output/) tab.

If you only want to run the variant calling steps, you can either manually select the variant calling processes or exclude downstream
processes (see: [Process selection](#process-selection)). Using manual process selection for example, you would run
`oncoanalyser` with the below command (assuming starting from FASTQ for DNA sequencing data):

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.neo_inputs.csv \
  --outdir output/ \
  --processes_manual \
  --processes_include alignment,redux,amber,cobalt,sage,pave,esvee,purple
```

### Why does `oncoanalyser` call too many / too few variants than another pipeline?

`oncoanalyser` uses variants with > 2% VAF. Other pipelines may have different assumptions which may cause differences in samples with low
tumor purity or a high number of subclonal variants.

### My HPC does not allow Docker

Some compute environments (especially HPCs) do not allow Docker as it runs a daemon as root which is deemed a security issue. In these
cases, using Singularity is recommended by providing `-profile singularity` when running `oncoanalyser` (also see:
[Setting up singularity](#setting-up-singularity)).

### Firewall prevents `oncoanalyser` from pulling Singularity images

Oncoanalyser can fail to pull Singularity images when there is a firewall set up. We recommend in these cases to manually download and cache
singularity containers. See [Downloading Apptainer containers](https://nf-co.re/docs/nf-core-tools/pipelines/download#downloading-apptainer-containers)
for more details.

### Network timeout

`oncoanalyser` may time out if pulling containers takes too long. To fix this, increase the network timeout in the config file
(see the [Nextflow config docs](https://www.nextflow.io/docs/latest/reference/config.html) to configure other container platforms):

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

Network timeout may also occur when downloading reference data. While the above solution might also work, we recommend downloading and
setting up reference data [manually](#manually-staging-reference-data) instead.

### Run fails due to insufficient CPUs/RAM/disk

Please see section: [Compute resources](#compute-resources)

### Can `oncoanalyser` CLI arguments be put in a config file?
Almost all `oncoanalyser` arguments in the [**Parameters tab**](https://nf-co.re/oncoanalyser/parameters) can be placed in a config file.

For example, the `oncoanalyser` arguments which start with `--` in this command:

```shell
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  -config refdata.config \
  --mode wgts \
  --genome GRCh38_hmf \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/outdir/
```

can be specified in a config file like so:
```groovy title='params.config'
params {
    mode = "wgts"
    genome = "GRCh38_hmf"
    input = "/path/to/samplesheet.csv"
    outdir = "/path/to/outdir/"
}
```

and provided as a config file when running `oncoanalyser`:
```shell
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.0.0 \
  -config refdata.config \
  -config params.config
```

:::tip

The `-config` Nextflow argument can be used multiple times to provide multiple config file

:::

### Errors and navigating the `work/` directory

When `oncoanalyser` crashes, you may need to further investigate error messages in the `.nextflow.log` files or the `work` directory.

The `work` directory contains the run scripts, logs, and input/output files for each process. Once the process is done running, only the
output files are 'published' (copied) to the final output directory (as specified by `--outdir`).

Below is an example `work` directory for one process. Error messages will typically be found in `.command.log`, `.command.err` or
`.command.out` log files. You can send these logs / error messages for example to the
`oncoanalyser` [Slack channel](https://nfcore.slack.com/channels/oncoanalyser), as an issue on the
[oncoanalyser](https://github.com/nf-core/oncoanalyser) or [WiGiTS](https://github.com/hartwigmedical/hmftools) GitHub.

```shell
work/
├── e5
│   └── f6e2e8f18ef70add9349164d5fb37e
│       ├── .command.sh     # Bash script used to run the process *within the container*
│       ├── .command.run    # Bash script used to run the process in the host machine
│       ├── .command.begin
│       ├── .command.log    # All log messages
│       ├── .command.err    # stderr log messages
│       ├── .command.out    # stdout log messages
│       ├── .command.trace  # Compute resource usage stats
│       ├── .exitcode       # Exit code
│       ├── <...>           # Input/output files or directories
│       └── versions.yml    # WiGiTS tool version
<...>
```

The `work/` directory can be hard to navigate due to the `<short_hash>/<long_hash>` structure. These hashes are shown (truncated) in the
console while running `oncoanalyser` (but can also be found in the `.nextflow.log` files):

```shell
[e5/f6e2e8] process > NFCORE_ONCOANALYSER:WGTS:REDUX_PROCESSING:REDUX (<group_id>_<sample_id>)     [100%] 2 of 2 ✔
```

Otherwise, you can use a utility like [tree](https://en.wikipedia.org/wiki/Tree_(command)) to show the directory structure, which allows
you to manually find the target process directory.

## Azure resource requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer to the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out
of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
