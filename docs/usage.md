# nf-core/oncoanalyser: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/oncoanalyser/usage](https://nf-co.re/oncoanalyser/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The `oncoanalyser` pipeline typically runs from FASTQ, BAM, or CRAM [input files](#analysis-starting-points), accepts
most GRCh37 and GRCh38 human [reference genome builds](#custom-genomes), and provides UMI ([unique molecular
identifier](#umi-processing)) processing for DNA sequencing data.

The pipeline supports two workflow modes: (1) whole genome and/or transcriptome, and (2) targeted panel. Both modes
accept DNA and RNA sequencing data from matched tumor / normal (with optional
[donor](#paired-tumor-and-normal-dna-with-donor-sample) sample) and tumor-only samples. The below table shows the
supported [sample setups](#sample-setups):

| Data Type | Tumor DNA          | Normal DNA         | Donor DNA          | Tumor RNA          |
| --------- | ------------------ | ------------------ | ------------------ | ------------------ |
| DNA       | :white_check_mark: | -                  | -                  | -                  |
| DNA       | :white_check_mark: | :white_check_mark: | -                  | -                  |
| DNA       | :white_check_mark: | :white_check_mark: | :white_check_mark: | -                  |
| DNA + RNA | :white_check_mark: | -                  | -                  | :white_check_mark: |
| DNA + RNA | :white_check_mark: | :white_check_mark: | -                  | :white_check_mark: |
| DNA + RNA | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| RNA       | -                  | -                  | -                  | :white_check_mark: |

## Running the pipeline

:::tip

Jump to [FAQ and troubleshooting](/oncoanalyser/2.1.0/docs/usage/faq_and_troubleshooting)

:::

A typical command for running `oncoanalyser` is shown below:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 2.1.0 \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

The [samplesheet](#samplesheet) provided to `--input` argument contains input sample details and corresponding files to
be analysed.

Additionally, various features of `oncoanalyser` can be configured by using a file provided to the `-config` argument.
This is generally recommended and it can be used to customise a number of settings or resources including:

- Reference genome and tool specific data: it is strongly recommended to [stage](#staging-reference-data) these files.
  Otherwise, `oncoanalyser` automatically stages them every run resulting in unnecessary disk/network usage
- Panel normalisation data: all panels except the built-in TSO500 panel require [additional
  setup](#panel-reference-data) of reference data
- [Other configuration](#custom-configuration): this may include [compute resources](#compute-resources) or [UMI
  settings](#umi-processing)

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
nextflow run nf-core/oncoanalyser -revision 2.1.0 -profile docker -params-file params.yaml
```

You can also generate such `yaml`/`json` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/oncoanalyser
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/oncoanalyser releases page](https://github.com/nf-core/oncoanalyser/releases) and find the latest pipeline version - numeric only (e.g. `2.1.0`). Then specify this when running the pipeline with `-revision` (one hyphen) - e.g. `-revision 2.1.0`. Of course, you can switch to another version by changing the number after the `-revision` flag.

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
| `filetype`      | File type: e.g. `fastq`, `bam`, `bai`; a full list of valid values can be found [here](https://github.com/nf-core/oncoanalyser/blob/2.1.0/lib/Constants.groovy#L56) |
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

#### BAM and CRAM

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

To run from CRAM, simply provide the CRAM and optionally the CRAM index with `bam` or `bai` in the `filetype` field:

```csv title="samplesheet.cram_crai.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.cram
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/other/dir/PATIENT1-T.dna.cram.crai
```

#### REDUX BAM

When running an analysis with DNA data from FASTQ, two of the most time consuming and resource intensive pipeline steps
are [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) read alignment and
[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) alignment processing. Where the REDUX output BAM
already exists for a given sample from a prior analysis, these read alignment and processing steps can be skipped by
providing the REDUX BAM as `bam_redux` in the `filetype` field. The REDUX BAM index can also optionally be provided with
`filetype` as `bai` if required.

```csv title="samplesheet.redux_bam_bai.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam_redux,/path/to/PATIENT1-T.dna.redux.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bai,/other/dir/PATIENT1-T.dna.redux.bam.bai
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

When starting from REDUX BAM, the filenames must have the format:

- `<sample_id>.redux.bam`
- `<sample_id>.redux.bam.bai`
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
and CLI arguments) to the Nextflow work directory during every run before proceeding with the analysis. It is therefore
strongly recommended to first stage and configure reference data to avoid repeated retrieval when performing multiple
`oncoanalyser` analyses.

#### Automatic staging

All reference data required for an analysis can be staged and prepared automatically by `oncoanalyser`. This is done by
configuring the desired analysis and then including the `--prepare_reference_only` argument, which causes `oncoanalyser`
to write reference data to the specified output directory without running the full pipeline.

For example the below samplesheet and command for analysing DNA data in `wgts` mode will stage the required `GRCh38_hmf`
genome (and indexes) and [WiGiTS](https://github.com/hartwigmedical/hmftools) resources files. As this analysis only
involves WGS data, no reference data files related to RNA or the `panel` mode will be retrieved.

```csv title="samplesheet.tn_dna.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.bam
```

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.1.0 \
  -profile docker \
  --mode wgts \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/ \
  --prepare_reference_only
```

Executing the above command will download and prepare default reference data without running any analysis, and once
complete the prepared reference files can be found in `./prepare_reference/reference_data/2.1.0/<datetimestamp>/`. You can then provide
a config file that points to these reference files (see [Configuring reference data](#configuring-reference-data)) which can
be used for subsequent `oncoanalyser` runs.

It is recommended to remove the Nextflow work directory once reference data staging is complete to free disk space.

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
  -revision 2.1.0 \
  -config refdata.config  \
  <...>
```

### Panel reference data

Analysis of panel / targeted sequencing data requires additional panel-specific reference data (e.g. region / gene
definitions, copy number and transcript normalisation data, known artefacts). This data is included and pre-configured
for the TSO500 panel, and can be used to analyse TSO500 sequence data by setting `--panel tso500` when running in
`targeted` mode:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.1.0 \
  -config refdata.config \
  -profile docker \
  --genome GRCh38_hmf \
  --mode targeted \
  --panel tso500 \
  --input samplesheet.csv \
  --outdir output/
```

For panels other than TSO500 (including whole exome), the panel-specific reference data must first be generated using a
training procedure detailed [here](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md).
The resulting panel-specific reference data must then be defined in a configuration file:

```groovy title="panel.config"
params {
    ref_data_panel_data_path = "/path/to/my_custom_panel_resources/"

    // These are relative paths within the dir provided by `ref_data_panel_data_path` above
    panel_data_paths {

        mycustompanel {  // This is the name that should be passed to the `--panel` argument

            // Genome version: '37' or '38'
            '38' {
                driver_gene_panel           = 'common/DriverGenePanel.custom_panel.38.tsv'
                sage_actionable_panel       = 'variants/ActionableCodingPanel.custom_panel.38.bed.gz'
                sage_coverage_panel         = 'variants/CoverageCodingPanel.custom_panel.38.bed.gz'
                pon_artefacts               = 'variants/pon_artefacts.custom_panel.38.tsv.gz'
                target_region_bed           = 'copy_number/target_regions_definition.custom_panel.38.bed.gz'
                target_region_normalisation = 'copy_number/cobalt_normalisation.custom_panel.38.tsv'
                target_region_ratios        = 'copy_number/target_regions_ratios.custom_panel.38.tsv'
                target_region_msi_indels    = 'copy_number/target_regions_msi_indels.custom_panel.38.tsv'

                // The below are optional and filepaths can be omitted for non-RNA panels by providing an empty list, e.g.:
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

To run an analysis of panel sequence data:

- provide both the panel-specific reference data configuration file via the `-config <file>` argument
- set the panel name in the `--panel <name>` argument, this must match the name defined in the configuration file
- set the `--force_panel` argument, which is required when not using the built-in `tso500` panel

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.1.0 \
  -config panel.config \
  -profile docker \
  --genome GRCh38_hmf \
  --mode targeted \
  --panel mycustompanel \
  --force_panel \
  --input samplesheet.csv \
  --outdir output/
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

Each index required for the analysis will first be created before running the rest of `oncoanalyser` with the following
command:

:::tip

In a process similar to [staging reference data](#automatic-staging), you can first generate the required indexes by
setting `--prepare_reference_only` and then provide the prepared reference files to `oncoanalyser` through a custom
config file. This avoids having to regenerate indexes for each new analysis.

:::

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.1.0 \
  -profile docker \
  -config genome.custom.config \
  --mode wgts \
  --genome CustomGenome \
  --genome_version <37|38> \
  --genome_type <alt|no_alt> \
  --force_genome \
  --input samplesheet.csv \
  --outdir output/
```

Creation of a STAR index also requires transcript annotations, please provide either of the following GTF files via the
`--ref_data_genome_gtf` option after decompressing:

- GRCh37: [GENCODE v37 (Ensembl v74)
  annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz)
- GRCh38: [GENCODE v38 (Ensembl v104)
  annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)

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

## Process selection

It is possible to exclude or include specific processes when running `oncoanalyser`. The full list of processes that can
be selected is available [here](https://github.com/nf-core/oncoanalyser/blob/2.1.0/lib/Constants.groovy#L32).

### Excluding processes

Most of the major components in `oncoanalyser` can be skipped using the `--processes_exclude` argument. There are
circumstances where it is desirable to skip resource intensive processes like VIRUSBreakend or where you have no use for
the outputs from some process such as the ORANGE report. In the example of skipping the VIRUSBreakend and ORANGE
processes, the `oncoanalyser` command would take the following form:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.1.0 \
  -profile docker \
  --mode wgts \
  --processes_exclude virusinterpreter,orange \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

:::warning

When skipping components no checks are done to identify orphan processes in the execution DAG or for redundant
processes.

:::

### Manual process selection

The `--processes_manual` argument can be used to enable manual process selection and `--processes_include
<process_1,process_2>` to configure individual processes to execute. One use case would be to run processes which are
not run by default, such as neoepitope calling with [NEO](https://github.com/hartwigmedical/hmftools/tree/master/neo).
To do this, provide the below example samplesheet:

```csv title='samplesheet.manual.csv'
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
PATIENT1,PATIENT1,PATIENT1-N,normal,dna,bam,/path/to/PATIENT1-N.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-T,tumor,dna,bam,/path/to/PATIENT1-T.dna.wgs.bam
PATIENT1,PATIENT1,PATIENT1-T-RNA,tumor,rna,bam,/path/to/PATIENT1-T.rna.wgs.bam
```

Then, run `oncoanalyser` with the `neo` process selected as well as all required upstream processes:

```bash
nextflow run nf-core/oncoanalyser \
  -revision 2.1.0 \
  -profile docker \
  --mode wgts \
  --processes_manual \
  --processes_include isofox,redux,amber,cobalt,sage,pave,esvee,purple,linx,lilac,neo \
  --genome GRCh38_hmf \
  --input samplesheet.neo_inputs.csv \
  --outdir output/
```

:::warning

It is the user's responsibility to select the required upstream processes for a downstream process to run. If not all
required processes are selected, `oncoanalyser` will not raise an error but instead finish without the downstream
process running.

:::

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
  -revision 2.1.0 \
  -profile docker \
  --mode wgts \
  --processes_manual \
  --processes_include neo \
  --genome GRCh38_hmf \
  --input samplesheet.neo_inputs.csv \
  --outdir output/
```

:::warning

Providing existing inputs will cause `oncoanalyser` to skip the corresponding process but none of the upstream
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
[nf-core documentation](https://nf-co.re/usage/configuration). Below subsections describe common use cases for custom configuration.

### Compute resources

The default compute resources (e.g. CPUs, RAM, disk space) configured in `oncoanalyser` may not be sufficient for one or
more processes. To change the resource requests, please see the [tuning workflow
resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) and [max
resources](https://nf-co.re/docs/usage/configuration#max-resources) sections of the nf-core website.

Below are the settings per WiGiTS tool that Hartwig uses internally and recommends. For high depth samples (e.g. panel
samples), you may need increase the memory for alignment, read processing (REDUX) and/or variant calling (SAGE or ESVEE)
steps.

```groovy
process {
    withName: '.*ALIGN'        { cpus = 12; memory = 72.GB; }
    withName: AMBER            { cpus = 16; memory = 24.GB; }
    withName: BAMTOOLS         { cpus = 16; memory = 24.GB; }
    withName: CHORD            { cpus = 4;  memory = 12.GB; }
    withName: COBALT           { cpus = 16; memory = 24.GB; }
    withName: CUPPA            { cpus = 4;  memory = 16.GB; }
    withName: 'ESVEE.*'        { cpus = 32; memory = 64.GB; }
    withName: LILAC            { cpus = 16; memory = 24.GB; }
    withName: 'LINX.*'         { cpus = 16; memory = 16.GB; }
    withName: REDUX            { cpus = 32; memory = 64.GB; }
    withName: ORANGE           { cpus = 4;  memory = 16.GB; }
    withName: 'PAVE.*'         { cpus = 8;  memory = 32.GB; }
    withName: PURPLE           { cpus = 8;  memory = 40.GB; }
    withName: 'SAGE.*'         { cpus = 32; memory = 64.GB; }
    withName: VIRUSBREAKEND    { cpus = 8;  memory = 64.GB; }
    withName: VIRUSINTERPRETER { cpus = 2;  memory = 8.GB;  }
}
```

Lastly, we recommend setting an upper limit on total resources that `oncoanalyser` is allowed to use. This will
typically be the max resources available to the VM / compute job. Below are the settings that Hartwig Medical Foundation
uses internally. When running multiple steps and/or samples in parallel, this will prevent `oncoanalyser` from
requesting more resources than available on the machine.

```groovy
process {
    resourceLimits = [
        cpus:   64,
        memory: 124.GB, // = 0.97 * 128.GB
        disk:   1500.GB,
        time:   48.h
    ]
}
```

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
