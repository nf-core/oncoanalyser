# nf-core/oncoanalyser: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/oncoanalyser/usage](https://nf-co.re/oncoanalyser/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The `oncoanalyser` pipeline typically runs from FASTQ, BAM, or CRAM [input files](#analysis-starting-points), accepts
most GRCh37 and GRCh38 human [reference genome builds](#custom-genomes), and provides UMI ([unique molecular
identifier](#umi-processing)) processing for DNA sequencing data.

Two main analysis modes are supported by `oncoanalyser`:

- [**wgts**](#whole-genome--transcriptome-sequencing-wgts): whole genome and/or transcriptome sequencing
- [**targeted**](#targeted-sequencing): targeted/panel sequencing

Both modes accept various combinations of DNA and/or RNA sequencing data from tumor-only or matched tumor / normal (with optional
[donor](#paired-tumor-and-normal-dna-with-donor-sample) sample). The below table shows the supported [sample setups](#sample-setups):

| DNA samples              | RNA samples |
| ------------------------ | ----------- |
| `tumor`                  | -           |
| `tumor`+`normal`         | -           |
| `tumor`+`normal`+`donor` | -           |
| `tumor`                  | `tumor`     |
| `tumor`+`normal`         | `tumor`     |
| `tumor`+`normal`+`donor` | `tumor`     |
| -                        | `tumor`     |

Besides the main analysis modes, several other modes are also available:

- [**purity_estimate**](#purity-estimate): tumor fraction estimation in longitudinal samples (e.g. for MRD)
- [**prepare_reference**](#automatic-staging): staging genomes and WiGiTS tool reference data
- [**panel_resource_creation**](#custom-panels): creating reference data for custom panels

## Running the pipeline

If you intend to run `oncoanalyser` on more than one sample, we recommend first [staging](#staging-reference-data) and
[configuring](#configuring-reference-data) reference data (genome and tool specific data). Otherwise, reference data is
automatically staged every run resulting in unnecessary disk/network usage.

A typical command for running `oncoanalyser` is shown below:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config reference_data.config \ # Optional but recommended
  -profile docker \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/ \
```

Below is a brief description of each argument:

- `-profile`: [configuration presets](#-profile) for different compute environments
- `-revision`: `oncoanalyser` version to run (can be a git [tag](https://github.com/nf-core/oncoanalyser/tags), [branch](https://github.com/nf-core/oncoanalyser/branches), or commit hash)
- `--mode`: [run mode](#run-modes)
- `--genome`: genome version, typically `GRCh38_hmf` or `GRCh37_hmf`
- `--input`: the [samplesheet](#samplesheet) containing sample details and corresponding files to be analysed
- `--output`: output directory
- `-config`: one or more configuration files for customising e.g. genome and tool specific data (as mentioned above),
  normalisation data for [custom panels](#custom-panels) (TSO500 panel supported by default), [compute resources](#compute-resources), or
  [other configuration](#custom-configuration)

:::tip

If you encounter any issues setting up or running `oncoanalyser`, please see
[FAQ and troubleshooting](/oncoanalyser/2.2.0/docs/usage/faq_and_troubleshooting)

:::

### Outputs

Running `oncoanalyser` will create the following files in your working directory:

```bash
work           # Directory containing the nextflow working files
<OUTDIR>       # Finished results in specified location (defined with --outdir)
.nextflow_log  # Log file from Nextflow
# Other nextflow hidden files, e.g. history of pipeline runs and old logs.
```

Descriptions of each output file in `<OUTDIR>` are provided in the [output](../output/) documentation.

### Reusing CLI arguments

To use the same CLI arguments across multiple runs, you can specify these in a `yaml` or `json` file via `-params-file <file>`.
The [above command](#running-the-pipeline) would have the equivalent `yaml` file:

```yaml title="params.yaml
mode: 'wgts'
genome: 'GRCh38_hmf'
input: 'samplesheet.csv'
outdir: 'output/'
<...>
```

and be run using this command:

```bash
nextflow run nf-core/oncoanalyser -revision 2.2.0 -profile docker -params-file params.yaml
```

You can also generate such `yaml`/`json` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/oncoanalyser
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/oncoanalyser releases page](https://github.com/nf-core/oncoanalyser/releases) and find the latest pipeline version - numeric only (e.g. `2.2.0`). Then specify this when running the pipeline with `-revision` (one hyphen) - e.g. `-revision 2.2.0`. Of course, you can switch to another version by changing the number after the `-revision` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, in the `<outdir>/pipeline_info/software_versions.yml` file.

To further assist in reproducibility, you can share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip

If you wish to share such a profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

:::

## Samplesheet

The samplesheet contains information in CSV format for each sample to be analysed by `oncoanalyser`, and uses a header
row as the first line with the below columns:

| Column          | Description                                                                                                                                                         |
| :-------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `group_id`      | Groups `sample_id` entries (e.g. tumor DNA, normal DNA, tumor RNA for one patient) into the same analysis                                                           |
| `subject_id`    | Subject/patient identifier, used internally to perform sanity check when processing multiple groups                                                                 |
| `sample_id`     | Sample identifier                                                                                                                                                   |
| `sample_type`   | Sample type: `tumor`, `normal`                                                                                                                                      |
| `sequence_type` | Sequence type: `dna`, `rna`                                                                                                                                         |
| `filetype`      | File type: e.g. `fastq`, `bam`, `bai`; a full list of valid values can be found [here](https://github.com/nf-core/oncoanalyser/blob/2.2.0/lib/Constants.groovy#L80) |
| `info`          | Additional sample information such as sequencing library and lane for [FASTQ](#fastq) files, this column is only required when running an analysis from FASTQ       |
| `filepath`      | Absolute filepath to input file, which can be a local filepath or supported protocol (http, https, ftp, s3, az, gz)                                                 |

The identifiers provided in the samplesheet are used to determine output file paths:

- `group_id`: top-level output directory for analysis files e.g. `output/PATIENT1/`
- tumor `sample_id`: output prefix for most filenames e.g. `PATIENT1-T.purple.sv.vcf.gz`
- normal `sample_id`: output prefix for some filenames e.g. `PATIENT1-N.cobalt.ratio.pcf`

### Analysis starting points

The `oncoanalyser` pipeline has been designed in such a way that allows an analysis to start from arbitrary entrypoints
as long as the required inputs are provided in the samplesheet. An analysis will generally start from either FASTQ or
alignment (BAM, CRAM, REDUX BAM) inputs, which are shown in the examples below.

#### FASTQ

To run from FASTQ:

- specify `fastq` in the `filetype` field,
- set sequencing library and lane information in the `info` field separated by `;`, and
- provide the forward ('R1') and reverse ('R2') FASTQ files in the `filepath` field separated by `;`

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,fastq,library_id:S1;lane:001,/path/to/PATIENT1-T_S1_L001_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L001_R2_001.fastq.gz
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,fastq,library_id:S1;lane:002,/path/to/PATIENT1-T_S1_L002_R1_001.fastq.gz;/path/to/PATIENT1-T_S1_L002_R2_001.fastq.gz
```

:::note

Currently only gzip compressed, non-interleaved paired-end FASTQ files are currently supported.

:::

#### BAM

To run from BAM, specify `bam` in the `filetype` field:

```csv title="samplesheet.bam.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
```

BAM indexes (.bai files) are expected to be in the same location as the BAM files with matching filenames suffixed with
`.bai`. Where this is not the case, you can also explicitly provide the BAM index location by specifying `bai` in the
`filetype` field:

```csv title="samplesheet.bam_bai.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/other/dir/PATIENT1-T.dna.bam.bai
```

#### CRAM

:::info

To run analyses starting from CRAM, you must use the CRAM format version ≤3.0 with the reference fully embedded. An
example command converting to the appropriate CRAM format is shown:

```bash
samtools view \
  --cram \
  --output-fmt-option version=3.0 \
  --output-fmt-option embed_ref=1 \
  --output-fmt-option reference=/path/to/reference.fasta \
  --output sample.cram \
  --threads 4 \
  --write-index \
  sample.bam
```

:::

To run from CRAM, use `cram` and `crai` in the `filetype` field. `crai` only needs to be provided if the CRAM index is
not in the same directory as the CRAM file:

```csv title="samplesheet.cram.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,cram,/path/to/PATIENT1-T.dna.cram
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,crai,/other/dir/PATIENT1-T.dna.cram.crai
```

Similarly, for REDUX CRAMs, provide `cram_redux` and optionally `crai`:

```csv title="samplesheet.redux_cram.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,cram_redux,/path/to/PATIENT1-T.dna.redux.cram
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,crai,/other/dir/PATIENT1-T.dna.cram.crai
```

:::warning

There is a fixed performance cost associated with reading CRAM files. This means the time it takes to read large CRAMs
vs BAMs is similar, whereas reading small CRAMs can take significantly longer (>10x) than reading small BAMs. If you
have small CRAMs (e.g. <10GB), it will be faster to decompress the CRAM into a BAM, and then run `oncoanalyser` with
this BAM.

This performance issue is due to how CRAM reading is implemented in
[htsjdk](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/CRAMFileReader.java) (which is
used throughout the WiGiTS tools). We plan to address this issue in future releases of `oncoanalyser`.

:::

#### REDUX BAM / CRAM

When running an analysis with DNA data from FASTQ, two of the most time consuming and resource intensive pipeline steps
are [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) read alignment and
[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) alignment processing.

`oncoanalyser` can be run starting from REDUX BAMs or CRAMs if they already exist from a prior analysis.

For REDUX BAMs, provide `bam_redux`/`cram_redux` in the `filetype` field, and optionally the BAM/CRAM index to `bai`/`crai` (only required
if indexes are not in the same directory as the BAM/CRAM):

```csv title="samplesheet.redux_bam_bai.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam_redux,/path/to/PATIENT1-T.dna.redux.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/other/dir/PATIENT1-T.dna.redux.bam.bai
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,cram_redux,/path/to/PATIENT2-T.dna.redux.cram
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,crai,/other/dir/PATIENT2-T.dna.redux.cram.crai
```

The `*.jitter_params.tsv` and `*.ms_table.tsv.gz` REDUX output files are expected to be in the same directory as the
REDUX BAM, and are required to run [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage). If these files
are located elsewhere, their paths can be explicitly provided by specifying `redux_jitter_tsv` and `redux_ms_tsv`in the
`filetype` field:

```csv title="samplesheet.redux_inputs.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam_redux,/path/to/PATIENT1-T.dna.redux.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,redux_jitter_tsv,/other/dir/PATIENT1-T.dna.jitter_params.tsv
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,redux_ms_tsv,/path/dir/PATIENT1-T.dna.ms_table.tsv.gz
```

:::tip

You can also [start from existing inputs](#starting-from-existing-inputs) other than from REDUX BAM

:::

:::warning

When starting from REDUX BAM/CRAM, the filenames must have the format:

- `<sample_id>.redux.bam` or `<sample_id>.redux.cram`
- `<sample_id>.redux.bam.bai` or `<sample_id>.redux.cram.crai`
- `<sample_id>.jitter_params.tsv`
- `<sample_id>.ms_table.tsv.gz`

For example, if `sample_id` is `PATIENT1-T`, the BAM filename must be `PATIENT1-T.redux.bam` and not e.g.
`PATIENT1.redux.bam`

:::

### Sample setups

Providing `sample_type` and `sequence_type` in different combinations allows `oncoanalyser` to run in different sample
setups. The below samplesheet examples use BAM files but different sample setups can also be specified for FASTQ or CRAM
files.

#### Paired tumor and normal DNA

```csv title="samplesheet.tn_dna.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
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

#### Paired tumor and normal DNA with tumor-only RNA

```csv title="samplesheet.wgts.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.bam
```

#### Paired tumor and normal DNA with donor sample

Including a donor sample in some types of analyses can be beneficial (e.g. bone marrow transplant) as this allows for
germline variant subtraction using both the patient's normal sample and the bone marrow donor's normal sample.

To include a donor sample in an analysis, specify `donor` in the `sample_type` field with a unique sample identifier in
the `sample_id` field:

```csv title="samplesheet.tn_with_donor.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.bam
PATIENT1,PATIENT1,PATIENT1-D,donor,dna,bam,/path/to/PATIENT1-D.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
```

### Multiple samples

To run with multiple samples, specify a different `group_id` and `subject_id` for each desired grouping:

```csv title="samplesheet.batch.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT2,PATIENT2,PATIENT2-N,normal,dna,bam,/path/to/PATIENT2-N.dna.bam
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,bam,/path/to/PATIENT2-T.dna.bam
```

## Reference data

The reference data used in `oncoanalyser` corresponds to and includes reference genomes and their indexes,
[WiGiTS](https://github.com/hartwigmedical/hmftools) resources files, and panel specific resource files. Descriptions
for each file can be found on the [WiGiTS resource file
documentation](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_RESOURCES.md) page.

### Staging reference data

By default `oncoanalyser` will download the required pre-configured reference data (based on the provided samplesheet
and CLI arguments) to the Nextflow work directory during every run before proceeding with the analysis.

However, strongly recommended to first stage and configure reference data to avoid repeated retrieval when
performing multiple `oncoanalyser` analyses. See the below for instructions.

#### Automatic staging

The reference data required for running `oncoanalyser` can be staged automatically using
`--mode prepare_reference` and specifying `--ref_data_types`. The below example command will stage the required
`GRCh38_hmf` genome (and indexes) and [WiGiTS](https://github.com/hartwigmedical/hmftools) resources files for WGS
analysis from BAM.

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -profile docker \
  --mode prepare_reference \
  --ref_data_types wgs \
  --genome GRCh38_hmf \
  --outdir output/
```

Once the above commands complete, the stated reference data can be found in `<outdir>/reference_data/2.2.0`. You will
then need to provide a config file that points to these reference files (see [Configuring reference data](#configuring-reference-data))
which can be used for subsequent `oncoanalyser` runs. The Nextflow work directory can also be removed to free up disk
space.

The below table shows the possible values for `--ref_data_types`. Note that multiple can be provided as comma separated
list, e.g. `--ref_data_types wgs,dna_alignment`

| Value                              | Description                                                                               | Combination of                                            |
| :--------------------------------- | :---------------------------------------------------------------------------------------- | :-------------------------------------------------------- |
| `wgs`                              | Ref data for WGS analysis from BAM                                                        | `fasta`, `fai`, `dict`, `img`, `hmftools`, `gridss_index` |
| `wts`                              | Ref data for WTS analysis from BAM                                                        | `fasta`, `fai`, `dict`, `img`, `hmftools`                 |
| `targeted`                         | Ref data for targeted analysis from BAM                                                   | `fasta`, `fai`, `dict`, `img`, `hmftools`, `panel`        |
| `bwamem2_index` or `dna_alignment` | BWA-MEM2 index. Required if aligning DNA FASTQs                                           |                                                           |
| `star_index` or `rna_alignment`    | STAR index. Required if aligning RNA FASTQs                                               |                                                           |
| `gridss_index`                     | GRIDSS index. Required if running Virusbreakend/Virusinterpreter                          |                                                           |
| `hmftools`                         | [WiGiTS](https://github.com/hartwigmedical/hmftools) resources files                      |                                                           |
| `panel`                            | Panel ref data. Only TSO500 currently supported. Please also specify arg `--panel tso500` |                                                           |

#### Manual staging

Where [automatic staging](#automatic-staging) cannot be used, reference data files can be downloaded manually from the
links provided in [Reference data URLs](#reference-data-urls). Any tarball should be extracted after download (using
`tar -xzvf <file>.tar.gz`).

To use locally staged reference data, see [Configuring reference data](#configuring-reference-data).

#### Configuring reference data

For `oncoanalyser` to use locally staged (or custom) reference data, the relevant settings can be defined in a
configuration file:

```groovy title="refdata.local.config"
params {
    genomes {
        GRCh38_hmf {
            fasta         = "/path/to/GRCh38_masked_exclusions_alts_hlas.fasta"
            fai           = "/path/to/GRCh38_masked_exclusions_alts_hlas.fasta.fai"
            dict          = "/path/to/GRCh38_masked_exclusions_alts_hlas.fasta.dict"
            img           = "/path/to/GRCh38_masked_exclusions_alts_hlas.fasta.img"
            bwamem2_index = "/path/to/bwa-mem2_index/"
            gridss_index  = "/path/to/gridss_index/"
            star_index    = "/path/to/star_index/"
        }
    }
    ref_data_hmf_data_path   = "/path/to/hmftools_data/"
    ref_data_panel_data_path = "/path/to/tso500_panel_data/"
}
```

The configuration file can then be supplied to `oncoanalyser` via the `-config <file>` argument:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config reference_data.config  \
  <...>
```

### Custom genomes

It is strongly recommended to use a Hartwig-distributed reference genome for alignments and subsequent analysis
(`GRCh37_hmf` or `GRCh38_hmf`). Where it is not feasible to do so, a custom genome can instead be used by providing the
relevant FASTA file in a configuration file:

:::warning

For GRCh38 genome builds, HLA typing and variant calling in `oncoanalyser` is incompatible with BAMs with fragments
aligned to HLA class I ALT contigs. These contigs should be removed or hard masked from the genome prior to use in
`oncoanalyser`. For cohorts with read data already mapped to a genome with HLA class I ALT contigs, alignments can
either be converted to FASTQ and provided to `oncoanalyser` with an appropriate genome build, or you can use
[Bamtools](https://github.com/hartwigmedical/hmftools/tree/master/bam-tools#altcontigremapper) to realign HLA reads to
the main assembly contigs.

:::

```groovy title='genome.custom.config'
params {
    genomes {
        CustomGenome {
            fasta = "/path/to/custom_genome.fa"
        }
    }
}
```

Each index can then be created in by using `--mode prepare_reference` and `--ref_data_types`
(see section [staging reference data](#automatic-staging)). The below example command would create the indexes for WGS analysis:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config genome.custom.config \
  -profile docker \
  --mode prepare_reference \
  --ref_data_types wgs,bwamem2_index,gridss_index
  --genome CustomGenome \
  --genome_version <37|38> \
  --genome_type <alt|no_alt> \
  --force_genome \
  --outdir output/
```

If aligning FASTQs from RNA seq data for WTS analysis, you should also provide `star_index` to `--ref_data_types`. Creating the
STAR index also requires transcript annotations; please provide either of the following GTF files via the `--ref_data_genome_gtf` option
after decompressing:

- GRCh37: [GENCODE v37 (Ensembl v74) annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz)
- GRCh38: [GENCODE v38 (Ensembl v104) annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)

:::warning

STAR index must use transcript annotations from Ensembl versions that match WiGiTS resource data (GRCh37: v74; GRCh38:
v104).

:::

When creating indexes for reference genomes with alternative haplotypes, an ALT file must be given with
`--ref_data_genome_alt`. Importantly, a STAR index will not be generated for reference genomes with alternative
haplotypes since this requires careful processing and is hence left to the user.

### Reference data URLs

_GRCh37 genome (Hartwig): `GRCh37_hmf`_

| Type                 | Link                                                                                                                                                                                                |
| :------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FASTA                | [Homo_sapiens.GRCh37.GATK.illumina.fasta](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/Homo_sapiens.GRCh37.GATK.illumina.fasta)                                      |
| FASTA index          | [Homo_sapiens.GRCh37.GATK.illumina.fasta.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai)          |
| FASTA seq dictionary | [Homo_sapiens.GRCh37.GATK.illumina.fasta.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict)        |
| BWA-MEM index image  | [Homo_sapiens.GRCh37.GATK.illumina.fasta.img](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/bwa_index_image/0.7.17-r1188/Homo_sapiens.GRCh37.GATK.illumina.fasta.img) |
| BWA-MEM2 index       | [bwa-mem2_index-2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/25.1/bwa-mem2_index-2.2.1.tar.gz)                                                              |
| GRIDSS index         | [gridss_index-2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/25.1/gridss_index-2.13.2.tar.gz)                                                                |
| STAR index           | [star_index-gencode_19-2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/25.1/star_index-gencode_19-2.7.3a.tar.gz)                                              |
| WiGiTS data          | [hmf_pipeline_resources.37_v2.1.0--1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.1.0--1.tar.gz)                            |
| TSO500 panel data    | [hmf_panel_resources.tso500.37_v2.0.0--3.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.37_v2.0.0--3.tar.gz)                      |

_GRCh38 genome (Hartwig): `GRCh38_hmf`_

| Type                 | Link                                                                                                                                                                                                  |
| :------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FASTA                | [GRCh38_masked_exclusions_alts_hlas.fasta](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/GRCh38_masked_exclusions_alts_hlas.fasta)                                      |
| FASTA index          | [GRCh38_masked_exclusions_alts_hlas.fasta.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/samtools_index-1.16/GRCh38_masked_exclusions_alts_hlas.fasta.fai)          |
| FASTA seq dictionary | [GRCh38_masked_exclusions_alts_hlas.fasta.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/samtools_index-1.16/GRCh38_masked_exclusions_alts_hlas.fasta.dict)        |
| BWA-MEM index image  | [GRCh38_masked_exclusions_alts_hlas.fasta.img](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/bwa_index_image-gatk-4.6.1.0/GRCh38_masked_exclusions_alts_hlas.fasta.img) |
| BWA-MEM2 index       | [bwa-mem2_index-2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/bwa-mem2_index-2.2.1.tar.gz)                                                                |
| GRIDSS index         | [gridss_index-2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/gridss_index-2.13.2.tar.gz)                                                                  |
| STAR index           | [star_index-gencode_38-2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/star_index-gencode_38-2.7.3a.tar.gz)                                                |
| WiGiTS data          | [hmf_pipeline_resources.38_v2.1.0--1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.1.0--1.tar.gz)                              |
| TSO500 panel data    | [hmf_panel_resources.tso500.38_v2.0.0--3.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.38_v2.0.0--3.tar.gz)                        |

## Run modes

### Whole genome / transcriptome sequencing (WGTS)

`--mode wgts` is used for analysing of whole genome (WGS) and/or whole transcriptome (WTS) sequencing data, and can be run like so:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config reference_data.config \
  -profile docker \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

### Targeted sequencing

`--mode targeted` is used for analysing targeted or panel sequencing samples. The TSO500 panel has in-built support by setting
`--panel tso500`. A typical run command for TSO500 panels would be:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config reference_data.config \
  -profile docker \
  --mode targeted \
  --panel tso500 \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

Panels other than TSO500 require additional arguments, as well as custom reference data to be created.
Please see [Custom panels](#custom-panels).

### Custom panels

`--mode panel_resource_creation` assists with creating custom panel reference data files (for panels other than TSO500), which fit and
normalise the biases inherent to that specific panel.

The below table summarises the required reference data files. Some panel reference data files must first be manually created - instructions
can be found on the [**WiGiTS targeted analysis readme**](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md).
Some these files are used with `--mode panel_resource_creation` to create the remaining required reference data files.

| Data type | File / config name            | Comment                                                                                                                 |
| :-------- | :---------------------------- | :---------------------------------------------------------------------------------------------------------------------- |
| DNA       | `driver_gene_panel`           | Manually created                                                                                                        |
| DNA       | `target_region_bed`           | Manually created                                                                                                        |
| DNA       | `target_region_msi_indels`    | Manually created                                                                                                        |
| DNA       | `target_region_ratios`        | Manually created                                                                                                        |
| DNA       | `target_region_normalisation` | Output from `--mode panel_resource_creation`                                                                            |
| DNA       | `pon_artefacts`               | Output from `--mode panel_resource_creation`                                                                            |
| RNA       | `isofox_gene_ids`             | Manually created                                                                                                        |
| RNA       | `isofox_tpm_norm`             | Output from `--mode panel_resource_creation`                                                                            |
| RNA       | `isofox_counts`               | Recommended to use `read_151_exp_counts.<ref_genome_version>.csv` from [WiGiTS reference data](#reference-data-urls)    |
| RNA       | `isofox_gc_ratios`            | Recommended to use `read_100_exp_gc_ratios.<ref_genome_version>.csv` from [WiGiTS reference data](#reference-data-urls) |

:::note

RNA reference data is only required if your panel supports RNA sequencing data.

:::

Once your manually created files are ready, create a samplesheet with a representative set of panel sequencing samples
(**≥20 recommended**). The below example samplesheet provides BAM files, but [FASTQ files](#fastq) can also be provided.

```csv title="samplesheet.panel_resource_creation.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/path/to/PATIENT1-T.dna.bam.bai
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,bam,/path/to/PATIENT2-T.dna.bam
PATIENT2,PATIENT2,PATIENT2-T,tumor,dna,bai,/path/to/PATIENT2-T.dna.bam.bai
```

Then, run `oncoanalyser` with `--mode panel_resource_creation` providing the samplesheet, as well as the relevant manually created files
to arguments `--driver_gene_panel`, `--target_regions_bed`, and `--isofox_gene_ids`:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config reference_data.config \
  -profile docker \
  --mode panel_resource_creation \
  --genome GRCh38_hmf \
  --input samplesheet.panel_resource_creation.csv \
  --driver_gene_panel DriverGenePanel.38.tsv \
  --target_regions_bed target_regions_definition.38.bed.gz \
  --isofox_gene_ids rna_gene_ids.csv \ # Optional, only provide if panel supports RNA sequencing data
  --outdir output/
```

Place the all the custom panel reference data files in a directory, and define the paths / file names in a configuration file:

```groovy title="panel.config"
params {
    ref_data_panel_data_path = "/directory/containing/my_custom_panel_resources/"

    // These are relative paths within the dir provided by `ref_data_panel_data_path` above
    panel_data_paths {

        my_custom_panel {  // This is the name that should be passed to the `--panel` argument

            // Genome version: '37' or '38'
            '38' {
                driver_gene_panel           = 'DriverGenePanel.38.tsv'
                pon_artefacts               = 'pave.somatic_artefacts.38.tsv'
                target_region_bed           = 'target_regions_definition.38.bed.gz'
                target_region_normalisation = 'cobalt.region_normalisation.38.tsv'
                target_region_ratios        = 'target_regions_ratios.38.tsv'
                target_region_msi_indels    = 'target_regions_msi_indels.38.tsv'

                // RNA. Optional, only provide if panel supports RNA data.
                isofox_gene_ids             = 'rna_gene_ids.csv'
                isofox_tpm_norm             = 'isofox.gene_normalisation.38.csv'
                isofox_counts               = 'read_151_exp_counts.37.csv'
                isofox_gc_ratios            = 'read_100_exp_gc_ratios.37.csv'
            }
        }
    }
}
```

Lastly, run `oncoanalyser` with `--mode targeted` to analyse your panel sequencing sample. You will also need to:

- provide the custom panel reference data configuration file to the `-config <file>` argument
- set the panel name in the `--panel <name>` argument as defined in the configuration file (e.g. `my_custom_panel`)
- set the `--force_panel` argument to enable non-built-in panels

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config reference_data.config \
  -config panel_data.config \
  -profile docker \
  --mode targeted \
  --panel my_custom_panel \
  --genome GRCh38_hmf \
  --force_panel \
  --input samplesheet.csv \
  --outdir output/
```

### Purity estimate

`--mode purity_estimate` uses [WISP](https://github.com/hartwigmedical/hmftools/tree/master/wisp) to estimate the tumor fraction
(aka purity) for a longitudinal sample (typically a ctDNA sample) guided by variants identified in a primary sample of the same patient
(typically a primary tissue biopsy). This can be used for example for detecting minimal residual disease (MRD).

The primary sample must first have been run in either [**WGTS**](#whole-genome--transcriptome-sequencing-wgts) or
[**targeted**](#targeted-sequencing) mode.

A samplesheet with the paths to the primary and longitudinal sample data is then created. Specifically:

- The BAM from the longitudinal tumor sample
- The AMBER and PURPLE directories from the **primary tumor** sample
- (Optional) The REDUX BAM of the normal sample, if the normal sample was provided in the primary sample run (i.e. was run in tumor/normal mode)

```csv title="samplesheet.purity_estimate.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
PATEINT1,PATIENT1,PATIENT1-L,tumor,dna,bam,longitudinal_sample,/path/to/PATIENT1-T.dna.longitudinal.bam
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam_redux,,/path/to/PATIENT1-N.dna.redux.bam
PATEINT1,PATIENT1,PATIENT1-T,tumor,dna,amber_dir,,/path/to/PATIENT1-T/amber/
PATEINT1,PATIENT1,PATIENT1-T,tumor,dna,purple_dir,,/path/to/PATIENT1-T/purple/
```

Then run `oncoanalyser` providing `--mode purity_estimate` and `--purity_estimate_mode <wgts|targeted>` (how the **longitudinal sample**
was sequenced):

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -config reference_data.config \
  -profile docker \
  --mode purity_estimate \
  --purity_estimate_mode targeted \
  --genome GRCh38_hmf \
  --input samplesheet.purity_estimate.csv \
  --outdir output/
```

:::note

`--purity_estimate_mode` simply sets different arguments for certain tools (e.g. SAGE). When running with `--purity_estimate_mode targeted`,
you do not need to configure panel ref data paths with as you would with `--mode targeted`.

:::

### Prepare reference data

`--mode prepare_reference` assists with staging all the reference data required to run `oncoanalyser`.
Please see: [Staging reference data: Automatic staging](#automatic-staging)

## Process selection

It is possible to exclude or manually select specific processes when running `oncoanalyser`. The full list of processes that can
be selected is available [here](https://github.com/nf-core/oncoanalyser/blob/2.2.0/lib/Constants.groovy#L53).

:::warning

It is the user's responsibility to select the required upstream processes for a downstream process to run. If not all
required processes are selected, `oncoanalyser` will not raise an error but instead finish without the downstream
process running.

:::

### Excluding processes

Most of the major components in `oncoanalyser` can be skipped using the `--processes_exclude` argument. You may want to
skip resource intensive processes like Virusbreakend, or ORANGE because you do not require the report, for example:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -profile docker \
  --mode wgts \
  --processes_exclude virusinterpreter,orange \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

### Manual process selection

The `--processes_manual` argument can be used to select the exact processes that `onconalyser` will run. For example,
you may only want to run alignment and SNV/indel, SV and CNV calling from DNA FASTQs, like so:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -profile docker \
  --mode wgts \
  --processes_manual alignment,redux,sage,amber,cobalt,esvee,sage,pave,purple \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

### Starting from existing inputs

An `oncoanalyser` analysis can start at arbitrary points as long as the required inputs are provided. For example,
neoepitope calling with [NEO](https://github.com/hartwigmedical/hmftools/tree/master/neo) can be run from existing
outputs generated by [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple),
[LILAC](https://github.com/hartwigmedical/hmftools/tree/master/lilac) and
[ISOFOX](https://github.com/hartwigmedical/hmftools/tree/master/isofox). To do this, provide the below example
samplesheet:

:::note

The original source input file (e.g. FASTQ, BAM, CRAM) must be provided for `oncoanalyser` to infer the correct analysis
type.

:::

```csv title='samplesheet.neo_inputs.csv'
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,isofox_dir,/path/to/PATIENT1.isofox_dir/
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,purple_dir,/path/to/PATIENT1.purple_dir/
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,linx_anno_dir,/path/to/PATIENT1.linx_anno_dir/
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,lilac_dir,/path/to/PATIENT1.lilac_dir/
```

Then, run `oncoanalyser` skipping all processes except for `neo`:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.2.0 \
  -profile docker \
  --mode wgts \
  --processes_manual neo \
  --genome GRCh38_hmf \
  --input samplesheet.neo_inputs.csv \
  --outdir output/
```

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
[nf-core documentation](https://nf-co.re/usage/configuration). Below subsections describe common use cases for custom configuration.

### Compute resources

The default compute resources (e.g. CPUs, RAM, disk space) configured in `oncoanalyser` may not be sufficient for one or
more processes (nf-core documentation: [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)).
For example, for high depth samples (e.g. panel samples), you may need increase the memory for alignment, read processing (REDUX),
small variant calling (SAGE), or structural variant calling (ESVEE) steps.

Below are the settings per tool that Hartwig Medical Foundation uses when running `oncoanalyser` in Google cloud:

```groovy
process {
    withName: '.*ALIGN'          { memory = 72.GB; cpus = 12; disk = 750.GB }
    withName: 'AMBER'            { memory = 24.GB; cpus = 16; disk = 375.GB }
    withName: 'BAMTOOLS'         { memory = 24.GB; cpus = 16; disk = 375.GB }
    withName: 'CHORD'            { memory = 12.GB; cpus = 4 ; disk = 375.GB }
    withName: 'CIDER'            { memory = 24.GB; cpus = 16; disk = 375.GB }
    withName: 'COBALT'           { memory = 24.GB; cpus = 16; disk = 375.GB }
    withName: 'CUPPA'            { memory = 16.GB; cpus = 4 ; disk = 375.GB }
    withName: 'ESVEE'            { memory = 96.GB; cpus = 32; disk = 375.GB }
    withName: 'ISOFOX'           { memory = 24.GB; cpus = 16; disk = 375.GB }
    withName: 'LILAC'            { memory = 24.GB; cpus = 16; disk = 375.GB }
    withName: 'LINX_.*'          { memory = 16.GB; cpus = 8 ; disk = 375.GB }
    withName: 'REDUX'            { memory = 64.GB; cpus = 32; disk = 750.GB }
    withName: 'ORANGE'           { memory = 16.GB; cpus = 4 ; disk = 375.GB }
    withName: 'PAVE.*'           { memory = 32.GB; cpus = 8 ; disk = 375.GB }
    withName: 'PEACH'            { memory = 4.GB ; cpus = 2 ; disk = 375.GB }
    withName: 'PURPLE'           { memory = 40.GB; cpus = 8 ; disk = 375.GB }
    withName: 'SAGE.*'           { memory = 64.GB; cpus = 32; disk = 375.GB }
    withName: 'TEAL.*'           { memory = 32.GB; cpus = 32; disk = 375.GB }
    withName: 'VIRUSBREAKEND'    { memory = 64.GB; cpus = 16; disk = 375.GB }
    withName: 'VIRUSINTERPRETER' { memory = 8.GB ; cpus = 2 ; disk = 375.GB }
    withName: 'WISP'             { memory = 16.GB; cpus = 4 ; disk = 375.GB }
}
```

We recommend setting an upper limit on total resources that `oncoanalyser` is allowed to use (nf-core
documentation: [max resources](https://nf-co.re/docs/usage/configuration#max-resources)). Otherwise, `oncoanalyser` may
crash when it tries to request more resources than available on a machine or compute job.
Below are some recommended resource limit settings:

```groovy
process {
    resourceLimits = [
        cpus:   64,
        memory: 120.GB, // Provides leeway on a 128.GB machine
        disk:   1500.GB,
        time:   48.h
    ]
}
```

The total runtime of `oncoanalyser` is ~3h for a paired 100x/30x tumor/normal WGS run starting from BAMs with parallel job execution via
Google batch. However, your runtime will vary depending on several factors such as sequencing depth, number of small/structural variants, or
parallel vs. non-parallel job execution.

### Container images

#### Custom containers

You may want to change which container or conda environment uses for a particular process (e.g. due to a newer tool
version being available). Please see [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) for
instructions.

#### Default containers

By default, `oncoanalyser` runs each tool using [Docker](https://docs.docker.com/engine/install/) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#) container images which are built by the
[Bioconda recipes](https://github.com/bioconda/bioconda-recipes/) CI/CD infrastructure. Below are links to these default
images should you want to download images manually (e.g. to [run `oncoanalyser`
offline](https://nf-co.re/docs/usage/getting_started/offline)).

**Docker (Bioconda)**

- Host: [quay.io](https://quay.io/organization/biocontainers)
- Repo URL example: https://quay.io/repository/biocontainers/hmftools-redux?tab=tags
- Image URI example: `quay.io/biocontainers/hmftools-redux:1.1--hdfd78af_1`

**Singularity (Bioconda)**

- Host: [Galaxy Project](https://depot.galaxyproject.org/singularity/)
- Image URI example: https://depot.galaxyproject.org/singularity/hmftools-redux:1.1--hdfd78af_1

**Bioconda recipes** for the above containers are found here:

- Host: [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes/)
- Recipe example: https://github.com/bioconda/bioconda-recipes/tree/master/recipes/hmftools-redux

**Docker** images built by Hartwig's CI/CD infrastructure are also available, intended for beta releases and not used by default in `oncoanalyser`

- Host: [Dockerhub](https://hub.docker.com/r/hartwigmedicalfoundation)
- Repo URL example: https://hub.docker.com/r/hartwigmedicalfoundation/redux/tags
- Image URI example: `docker.io/hartwigmedicalfoundation/redux:1.1`

:::tip

You can get the URIs for the default container images from the `oncoanalyser` repo with the below shell commands:

- Docker: `grep -rohE "'biocontainers.*'" oncoanalyser/modules/local/ | sort | uniq`
- Singularity: `grep -rohE "'https://depot.galaxyproject.*'" oncoanalyser/modules/local/ | sort | uniq`

:::

#### Container configuration

All configuration options for containers can be found in the [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/reference/config.html).
A typical config might look like this:

```groovy
singularity {
    enabled = true
    cacheDir = '/path/to/cache_dir/'
    autoMounts = true
    runOptions = "-B </path/to/desired/mounted/volume/>"
    pullTimeout = '2h'
}
```

### Executors

The [executor](https://www.nextflow.io/docs/latest/executor.html) is a Nextflow component that allows to submission of jobs for example via
[SLURM](https://www.nextflow.io/docs/latest/executor.html#slurm) (typically on an HPC),
[AWS Batch](https://www.nextflow.io/docs/latest/aws.html), or
[Google Batch](https://www.nextflow.io/docs/latest/google.html).

To enable SLURM for example, you would provide the below config:

```groovy
process {
    executor = "slurm"
}
```

Additional options for the enabled executor can be provided to the `executor` directive as shown below. See the
[Config: Executor](https://www.nextflow.io/docs/latest/reference/config.html#executor) Nextflow documentation for all options.

```groovy
executor {
    queueSize         = 100
    queueStatInterval = '10 sec'
    pollInterval      = '10 sec'
    submitRateLimit   = '10 sec'
}
```

### Custom tool arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines
provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the
[customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### UMI processing

Unique molecular identifiers (UMI) allow for read deduplication and error correction. UMI processing is performed by
[fastp](https://github.com/OpenGene/fastp?tab=readme-ov-file#unique-molecular-identifier-umi-processing) for FASTQ
files, and [REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux#deduplication) for BAM files. Depending
on the presence/format of your UMI strings, you may need to configure one or more of these arguments:

```groovy title='umi.config'
params {
    // For FASTQ files
    fastp_umi_enabled = true        // Enable UMI processing by fastp
    fastp_umi_location = "per_read" // --umi_loc fastp arg
    fastp_umi_length = 7            // --umi_len fastp arg
    fastp_umi_skip = 0              // --umi_skip fastp arg

    // For BAM files
    redux_umi_enabled = true        // Enable UMI processing by REDUX
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

## Azure resource requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer to the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out
of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted via your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
