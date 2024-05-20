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
- granular control over reference/resource files

These features enable oncoanalyser to be run in a highly flexible way. For example, an analysis can be run with existing
PURPLE data as the starting point and skip variant calling processes. Additionally, reference/resource files can be
staged locally to optimise execution or modified to create user-defined driver gene panels.

:::danger

When starting from BAMs rather than FASTQ it is expected that:

- RNA read alignments are generated with STAR using [specific
  parameters](https://github.com/hartwigmedical/hmftools/tree/master/isofox#a-note-on-alignment-and-multi-mapping), this
  is **critical** for WTS data, and
- reads are aligned to a Hartwig-distributed reference genome ([custom genomes](#custom-genomes) can be used but are not
  recommended)

:::

## Supported analyses

A variety of analyses are accessible in oncoanalyser and are implicitly run according to the data described in the
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

Several different input filetypes beyond FASTQ and BAM are recognised, including intermediate output files generated
during execution such as the PURPLE output directory. The full list of recognised input filetypes is available
[here](https://github.com/nf-core/oncoanalyser/blob/0.4.5/lib/Constants.groovy#L58-L90).

### Simple example

#### FASTQ

:::note

Currently only non-interleaved paired-end reads are accepted as FASTQ input.

:::

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
P1__wgts,P1,SA,normal,dna,fastq,library_id:SA_library;lane:001,/path/to/P1.SA.normal.dna.wgs.001.R1.fastq.gz;/path/to/P1.SA.normal.dna.wgs.001.R2.fastq.gz
P1__wgts,P1,SB,tumor,dna,fastq,library_id:SB_library;lane:001,/path/to/P1.SB.tumor.dna.wgs.001.R1.fastq.gz;/path/to/P1.SB.tumor.dna.wgs.001.R2.fastq.gz
P1__wgts,P1,SC,tumor,rna,fastq,library_id:SC_library;lane:001,/path/to/P1.SC.tumor.rna.wts.001.R1.fastq.gz;/path/to/P1.SC.tumor.rna.wts.001.R2.fastq.gz
```

#### BAM

:::note

Inputs with the `bam` filetype will be processed by MarkDups as required by hmftools. Where an input BAM has already
been processed specifically by [HMF
MarkDups](https://github.com/hartwigmedical/hmftools/blob/master/mark-dups/README.md), you can avoid needless
reprocessing by setting `bam_markdups` as the filetype instead. It is important to understand that duplicate marking by
other tools (e.g. GATK) cannot be used as a substitute since HMF MarkDups performs key operations beyond just duplicate
marking.

<br />

Please note there are other essential requirements around the use of BAMs as inputs, see the warning above in the
[Introduction](#introduction).

:::

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

| Column        | Description                                                                    |
| ------------- | ------------------------------------------------------------------------------ |
| group_id      | Group ID for a set of samples and inputs                                       |
| subject_id    | Subject/patient ID                                                             |
| sample_id     | Sample ID                                                                      |
| sample_type   | Sample type: `tumor`, `normal`                                                 |
| sequence_type | Sequence type: `dna`, `rna`                                                    |
| filetype      | File type: e.g. `fastq`, `bam`, `bai`                                          |
| info          | Additional input information: `library_id`, `lane`, `cancer_type` _[optional]_ |
| filepath      | Absolute filepath to input file (can be local filepath, URL, S3 URI)           |

The identifiers provided in the samplesheet are used to set output file paths:

- `group_id`: top-level output directory for analysis files e.g. `output/COLO829_example/`
- tumor `sample_id`: output prefix for most filenames e.g. `COLO829T.purple.sv.vcf.gz`
- normal `sample_id`: output prefix for some filenames e.g. `COLO829R.cobalt.ratio.pcf`

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 0.4.5 \
  --mode <wgts|targeted> \
  --genome <GRCh37_hmf|GRCh38_hmf> \
  --input samplesheet.csv \
  --outdir <output_directory>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information on profiles.

:::note

Reference data will be retrieved by oncoanalyser for every analysis run. It is therefore strongly recommended when
running multiple analyses to pre-stage reference data locally to avoid it being retrieved multiple times. See [Staging
reference data](#staging-reference-data).

:::

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

Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must
only be used for [tuning process resource
specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such
as output directories), or module arguments (args).

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

To further assist in reproducbility, you can share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Advanced usage

### Selecting processes

Most of the major components in oncoanalyser can be skipped using `--processes_exclude` (the full list of available
processes can be view [here](https://github.com/nf-core/oncoanalyser/blob/0.4.5/lib/Constants.groovy#L36-L56)).
Multiple processes can be given as a comma-separated list. While there are some use-cases for this feature (e.g.
skipping resource intensive processes such as VIRUSBreakend), it becomes more powerful when combined with existing
inputs as described in the following section.

:::warning

When skipping components no checks are done to identify orphan processes in the execution DAG or for redundant
processes.

:::

### Existing inputs

The oncoanalyser pipeline has been designed to allow entry at arbitrary points, which is particularly useful in
situations where previous outputs exist and re-running oncoanalyser is desired (e.g. to subsequently execute an
optional sensor or use an upgrade component such as PURPLE). The primary advantage of this approach is that only the
required processes are executed, reducing costs and runtimes by skipping unnecessary processes.

In order to effectively utilise this feature, existing inputs must be set in the [samplesheet](#samplesheet) and the
appropriate [processes selected](#selecting-processes). Take the below example where existing PURPLE inputs are used so
that all upstream variant calling can be skipped:

```csv title='samplesheet.existing_purple.csv'
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
P1__wgts,P1,SA,normal,dna,bam,/path/to/P1.SA.normal.dna.wgs.bam
P1__wgts,P1,SB,tumor,dna,bam,/path/to/P1.SB.tumor.dna.wgs.bam
P1__wgts,P1,SB,tumor,dna,purple_dir,/path/to/P1.purple_dir/
```

:::note

The original source input file (i.e. BAM or FASTQ) must always be provided for oncoanalyser to infer the correct
analysis type.

:::

And now run and skip variant calling:

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 0.4.5 \
  --mode wgts \
  --processes_exclude markdups,amber,cobalt,gridss,gripss,sage,pave \
  --genome GRCh38_hmf \
  --input samplesheet.csv \
  --outdir output/
```

:::warning

Providing existing inputs will cause oncoanalyser to skip the corresponding process but _not any_ of the upstream
processes. It is the responsibility of the user to skip all relevant processes.

:::

### Configuring reference data

All reference data can be configured as needed, and are defined in the following locations:

| Reference data          | Filepath                  | Note                                    |
| ----------------------- | ------------------------- | --------------------------------------- |
| Genomes and indexes     | `conf/hmf_genomes.config` | Absolute paths                          |
| hmftools resource files | `conf/hmf_data.config`    | Paths relative to data bundle directory |
| Panel resource files    | `conf/panel_data.config`  | Paths relative to data bundle directory |

See the below sections for further details on customising reference data.

#### Customising hmf data

To override hmftools resource files, first [stage the bundle](#staging-reference-data) locally then copy in your
custom file under the bundle directory and create a new config with relevant file paths:

```text title="hmf_data.custom.config"
params {
    hmf_data_paths {
        '38' {
            driver_gene_panel     = 'custom_files/DriverGenePanel.tsv'
            sage_actionable_panel = 'custom_files/ActionableCodingPanel.bed.gz'
            sage_coverage_panel   = 'custom_files/CoverageCodingPanel.bed.gz'
        }
    }
}
```

To use these hmftools resource file overrides in oncoanalyser the local bundle directory must be provided with
`--ref_data_hmf_data_path`.

#### Customise other data

The path or URI to the VIRUSBreakend database can also be explicitly set with `--ref_data_virusbreakenddb_path`. There
are additional arguments to manually set various other reference data files, please review the parameters documentation
for the complete list.

#### Staging reference data

Default reference data can be staged locally with oncoanalyser by providing a samplesheet for the desired analysis and
setting the `--prepare_reference_only` argument. The samplesheet and oncoanalyser configuration will determine the
relevant reference data to download. For example the following command will download the `GRCh38_hmf` genome plus
indices, reference data, and databases required to run a WGTS analysis for tumor/normal DNA with tumor RNA:

```csv title="samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
P1__wgts,P1,SA,normal,dna,bam,/path/to/P1.SA.normal.dna.wgs.bam
P1__wgts,P1,SB,tumor,dna,bam,/path/to/P1.SB.tumor.dna.wgs.bam
P1__wgts,P1,SC,tumor,rna,bam,/path/to/P1.SC.tumor.rna.wts.bam
```

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 0.4.5 \
  --mode wgts \
  --genome GRCh38_hmf \
  --prepare_reference_only \
  --input samplesheet.csv \
  --outdir prepare_reference/
```

Executing the above command will download and unpack default reference data without running any analysis, and once
complete the prepared reference files can found in `./prepare_reference/reference_data/0.4.5/<datetimestamp>/`. It is
recommended to remove the Nextflow work directory after staging data to free disk space.

For oncoanalyser to use locally staged reference data a custom config can be used:

```text title="refdata.local.config"
params {

    genomes {
        GRCh38_hmf {
            fasta           = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
            fai             = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
            dict            = "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict"
            bwamem2_index   = "/path/to/bwa-mem2_index/"
            gridss_index    = "/path/to/gridss_index/"
            star_index      = "/path/to/star_index/"
        }
    }

    ref_data_hmf_data_path        = "/path/to/hmftools_data/"
    ref_data_panel_data_path      = "/path/to/tso500_panel_data/"
    ref_data_virusbreakenddb_path = "/path/to/virusbreakenddb/"
}
```

Specific reference files can also be downloaded directly from the hosting service with the corresponding URL.

##### Reference data URLs

_GRCh37 genome (Hartwig) [`GRCh37_hmf`]_

| Type                 | Name                                                                                                                                                                                         |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FASTA                | [Homo_sapiens.GRCh37.GATK.illumina.fasta](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/Homo_sapiens.GRCh37.GATK.illumina.fasta)                               |
| FASTA index          | [Homo_sapiens.GRCh37.GATK.illumina.fasta.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai)   |
| FASTA seq dictionary | [Homo_sapiens.GRCh37.GATK.illumina.fasta.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict) |
| bwa-mem2 index       | [bwa-mem2_index/2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                       |
| GRIDSS index         | [gridss_index/2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                         |
| STAR index           | [star_index/gencode_19/2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/star_index/gencode_19/2.7.3a.tar.gz)                                       |

_GRCh38 genome (Hartwig) [`GRCh38_hmf`]_

| Type                 | Name                                                                                                                                                                                                         |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| FASTA                | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna)                               |
| FASTA index          | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai)   |
| FASTA seq dictionary | [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/samtools_index/1.16/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict) |
| bwa-mem2 index       | [bwa-mem2_index/2.2.1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/bwa-mem2_index/2.2.1.tar.gz)                                                                       |
| GRIDSS index         | [gridss_index/2.13.2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.1/gridss_index/2.13.2.tar.gz)                                                                         |
| STAR index           | [star_index/gencode_38/2.7.3a.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/24.0/star_index/gencode_38/2.7.3a.tar.gz)                                                       |

_Other reference data_

| Type                   | Name                                                                                                                                           |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| hmftools data (GRCh37) | [hmftools/5.34_37--2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/5.34_37--2.tar.gz)                |
| hmftools data (GRCh38) | [hmftools/5.34_38--2.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/5.34_38--2.tar.gz)                |
| TSO500 data (GRCh37)   | [panels/tso500_5.34_37--1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/tso500_5.34_37--1.tar.gz)      |
| TSO500 data (GRCh38)   | [panels/tso500_5.34_38--1.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/tso500_5.34_38--1.tar.gz)      |
| HLA slice BED          | [hla_slice/grch38_alt.plus_homologous.bed](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/other/hla_slice/grch38_alt.plus_homologous.bed) |
| VIRUSBreakend database | [virusbreakenddb_20210401.tar.gz](https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/virusbreakend/virusbreakenddb_20210401.tar.gz)           |

#### Custom genomes

It is strongly recommended to use a Hartwig-distributed reference genome for alignments and subsequent analysis
(`GRCh37_hmf` or `GRCh38_hmf`). Where it is not feasible to do so, a custom genome can instead be used by providing the
relevant FASTA file in a configuration file:

```text title='genome.custom.config'
params {
    genomes {
        CustomGenome {
            fasta = "/path/to/custom_genome.fa"
        }
    }
}
```

Each index required for the analysis will first be created before running the rest of oncoanalyser with the following
command:

:::note

In a process similar to [staging reference data](#staging-reference-data), you can first generate the required indexes
by setting `--prepare_reference_only` and then provide the prepared reference files to oncoanalyser through a custom
config file. This avoids having to regenerate indexes for each new analysis.

:::

```bash
nextflow run nf-core/oncoanalyser \
  -profile docker \
  -revision 0.4.5 \
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

Creation of a STAR index also requires transcript annotations, please provide either of the following GTF files via the
`--ref_data_genome_gtf` option after decompressing:

- GRCh37: [GENCODE v38 (Ensembl v104)
  annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)
- GRCh38: [GENCODE v37 (Ensembl v74)
  annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz)

:::warning

STAR index must use transcript annotations from Ensembl versions that match hmftools resource data (GRCh37: v74; GRCh38:
v104).

:::

When creating indexes for reference genomes with alternative haplotypes, an ALT file must be given with
`--ref_data_genome_alt`. Importantly, a STAR index will not be generated for reference genomes with alternative
haplotypes since this requires careful processing and is hence left to the user.

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

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

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

### Custom containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom tool arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure resource requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer to the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

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
