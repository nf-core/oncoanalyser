# nf-core/oncoanalyser: Output

## Introduction

This document describes the output produced by `oncoanalyser`. The pipeline writes output files to the below directory tree structure.
Files are grouped by the `group_id` provided in the sample sheet, then by tool.

```shell
output/
│
├── group_id_1/
│   ├── alignments/
│   ├── amber/
│   ├── bamtools/
│   ├── chord/
│   ├── cobalt/
│   ├── cuppa/
│   ├── esvee/
│   ├── isofox/
│   ├── lilac/
│   ├── linx/
│   ├── orange/
│   ├── pave/
│   ├── purple/
│   ├── sage/
│   ├── sigs/
│   ├── virusbreakend/
│   └── virusinterpreter/
│
├── group_id_2/
│   └── ...
│
...
│
└── pipeline_info/
```

## Pipeline overview
- [Read alignment](#read-alignment)
  - [BWA-MEM2](#bwa-mem2) - DNA read alignment
  - [STAR](#star) - RNA read alignment
- [Read post-processing](#read-post-processing)
  - [REDUX](#redux) - DNA read post-processing
  - [Picard MarkDuplicates](#picard-markduplicates) - RNA read post-processing
- [SNV, MNV, INDEL calling](#snv-mnv-indel-calling)
  - [SAGE](#sage) - Small variant calling
  - [PAVE](#pave) - Transcript/coding effect annotation
- [SV calling](#sv-calling)
  - [ESVEE](#esvee) - Structural variant calling
- [CNV calling](#cnv-calling)
  - [AMBER](#amber) - B-allele frequencies
  - [COBALT](#cobalt) - Read depth ratios
  - [PURPLE](#purple) - Purity/ploidy estimation, small variant annotation
- [SV and driver event interpretation](#sv-and-driver-event-interpretation)
  - [LINX](#linx) - SV and driver event interpretation
- [RNA transcript analysis](#rna-transcript-analysis)
  - [ISOFOX](#isofox) - RNA transcript analysis
- [Oncoviral detection](#oncoviral-detection)
  - [VIRUSBreakend](#virusbreakend) - Viral content and integration calling
  - [VirusInterpreter](#virusinterpreter) - Post-processing
- [Immune analysis](#immune-analysis)
  - [LILAC](#lilac) - HLA typing
  - [NEO](#neo) - Neo-epitope prediction
- [Mutational signature fitting](#mutational-signature-fitting)
  - [SIGS](#sigs) - Mutational signature fitting
- [HRD prediction](#hrd-prediction)
  - [CHORD](#chord) - HRD prediction
- [Tissue of origin prediction](#tissue-of-origin-prediction)
  - [CUPPA](#cuppa) - Tissue of origin prediction
- [Summary report](#summary-report)
  - [ORANGE](#orange) - Summary report
- [Pipeline information](#pipeline-information) - Pipeline execution metrics

## Pipeline outputs

### Read alignment

#### BWA-MEM2

[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) is a short-read mapping tool used to align reads to a large reference sequences.
In `Oncoanalyser`, it is used to align DNA reads to the human genome.

_No outputs are published directly from bwa-mem2, see [REDUX](#redux) for the fully processed alignment outputs_

#### STAR

[STAR](https://github.com/alexdobin/STAR) is a specialised mapping used to align RNA reads to a reference transcriptome.

_No outputs are published directly from STAR, see [Picard MarkDuplicates](#picard-markduplicates) for the fully processed alignment outputs_

### Read post-processing

#### REDUX

[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) applies various alignment post-processing routines such as duplicate
marking and unmapping of problematic regions for DNA reads. It can also handle UMIs when configured to do so.

```shell
<group_id>/alignments/
├── dna
│   ├── <tumor_dna_id>.jitter_params.tsv          # Tumor DNA sample: Microsatellite jitter model parameters
│   ├── <tumor_dna_id>.ms_table.tsv.gz            # Tumor DNA sample: Aggregated repeat units and repeat counts
│   ├── <tumor_dna_id>.redux.bam                  # Tumor DNA sample: Read alignments
│   ├── <tumor_dna_id>.redux.bam.bai              # Tumor DNA sample: Read alignments index
│   ├── <tumor_dna_id>.redux.duplicate_freq.tsv   # Tumor DNA sample: Duplicate read frequencies
│   ├── <tumor_dna_id>.repeat.tsv.gz              # Tumor DNA sample: Repeat units and repeat counts per site
│   ├── <normal_dna_id>.jitter_params.tsv         # Normal DNA sample: Microsatellite jitter model parameters
│   ├── <normal_dna_id>.ms_table.tsv.gz           # Normal DNA sample: Aggregated repeat units and repeat counts
│   ├── <normal_dna_id>.redux.bam                 # Normal DNA sample: Read alignments
│   ├── <normal_dna_id>.redux.bam.bai             # Normal DNA sample: Read alignments index
│   ├── <normal_dna_id>.redux.duplicate_freq.tsv  # Normal DNA sample: Duplicate read frequencies
│   └── <normal_dna_id>.repeat.tsv.gz             # Normal DNA sample: Repeat units and repeat counts per site
```

#### Picard MarkDuplicates

[Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) used to mark duplicate reads
following alignment of RNA reads in `oncoanalyser`.

```shell
└── rna
    ├── <tumor_rna_id>.md.bam                     # Tumor DNA sample: Read alignments
    ├── <tumor_rna_id>.md.bam.bai                 # Tumor DNA sample: Read alignments index
    └── <tumor_rna_id>.md.metrics                 # Tumor DNA sample: Duplicate read marking metrics
```

### SNV, MNV, INDEL calling

#### SAGE

[SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage) is an SNV, MNV, and INDEL caller optimised for 100x tumor and 40x normal.

```shell
<group_id>/sage
├── somatic
│   ├── <normal_dna_id>.sage.bqr.png             # Normal DNA sample: Base quality recalibration metrics plot
│   ├── <normal_dna_id>.sage.bqr.tsv             # Normal DNA sample: Base quality recalibration metrics
│   ├── <tumor_dna_id>.sage.bqr.png              # Tumor DNA sample: Base quality recalibration metrics plot
│   ├── <tumor_dna_id>.sage.bqr.tsv              # Tumor DNA sample: Base quality recalibration metrics
│   ├── <tumor_dna_id>.sage.exon.medians.tsv     # Tumor DNA sample: Exon median depths
│   ├── <tumor_dna_id>.sage.gene.coverage.tsv    # Tumor DNA sample: Gene coverages
│   ├── <tumor_dna_id>.sage.somatic.vcf.gz       # Tumor DNA sample: Small variant calls
│   └── <tumor_dna_id>.sage.somatic.vcf.gz.tbi   # Tumor DNA sample: Small variant calls index
├── germline
│   ├── <normal_dna_id>.sage.bqr.png             # Normal DNA sample: Base quality recalibration metrics plot
│   ├── <normal_dna_id>.sage.bqr.tsv             # Normal DNA sample: Base quality recalibration metrics
│   ├── <normal_dna_id>.sage.exon.medians.tsv    # Normal DNA sample: Exon median depths
│   ├── <normal_dna_id>.sage.gene.coverage.tsv   # Normal DNA sample: Gene coverages
│   ├── <tumor_dna_id>.sage.bqr.png              # Tumor DNA sample: Base quality recalibration metrics plot
│   ├── <tumor_dna_id>.sage.bqr.tsv              # Tumor DNA sample: Base quality recalibration metrics
│   ├── <tumor_dna_id>.sage.germline.vcf.gz      # Tumor DNA sample: Germline small variant calls
│   └── <tumor_dna_id>.sage.germline.vcf.gz.tbi  # Tumor DNA sample: Germline small variant calls index
└── append
    ├── <normal_dna_id>.sage.append.vcf.gz       # Normal DNA sample: VCF with SMNVs and RNA data appended
    └── <tumor_dna_id>.sage.append.vcf.gz        # Tumor DNA sample: VCF with SMNVs and RNA data appended
```

#### PAVE

[PAVE](https://github.com/hartwigmedical/hmftools/tree/master/pave) annotates variants called by SAGE with impact information with regards
to transcript and coding effects.

```shell
<group_id>/pave/
├── <tumor_dna_id>.sage.germline.pave.vcf.gz      # Tumor DNA sample: VCF with annotated germline SMNVs
├── <tumor_dna_id>.sage.germline.pave.vcf.gz.tbi  # Tumor DNA sample: VCF index
├── <tumor_dna_id>.sage.somatic.pave.vcf.gz       # Tumor DNA sample: VCF with annotated somatic SMNVs
└── <tumor_dna_id>.sage.somatic.pave.vcf.gz.tbi   # Tumor DNA sample: VCF index
```

### SV calling

#### ESVEE

[ESVEE](https://github.com/hartwigmedical/hmftools/tree/master/esvee) is a structural variant caller than uses both read support and local
breakend/breakpoint assemblies to call variants.

```shell
<group_id>/esvee/
├── prep
│   ├── <tumor_dna_id>.esvee.prep.bam                  # Tumor DNA sample: BAM with candidate SV reads
│   ├── <tumor_dna_id>.esvee.prep.bam.bai              # Tumor DNA sample: BAM index
│   ├── <tumor_dna_id>.esvee.prep.disc_stats.tsv       # Tumor DNA sample: Discordant reads stats
│   ├── <tumor_dna_id>.esvee.prep.fragment_length.tsv  # Tumor DNA sample: Fragment length stats
│   ├── <tumor_dna_id>.esvee.prep.junction.tsv         # Tumor DNA sample: Candidate junctions
│   ├── <normal_dna_id>.esvee.prep.bam                 # Tumor DNA sample: BAM with candidate SV reads
│   └── <normal_dna_id>.esvee.prep.bam.bai             # Tumor DNA sample: BAM index
├── assemble
│   ├── <tumor_dna_id>.esvee.assembly.tsv              # Tumor DNA sample: Breakend assemblies
│   ├── <tumor_dna_id>.esvee.alignment.tsv             # Tumor DNA sample: Assemblies realigned to the ref genome
│   ├── <tumor_dna_id>.esvee.breakend.tsv              # Tumor DNA sample: Breakends
│   ├── <tumor_dna_id>.esvee.phased_assembly.tsv       # Tumor DNA sample: Phased assemblies
│   ├── <tumor_dna_id>.esvee.raw.vcf.gz                # Tumor DNA sample: VCF with candidate breakends
│   └── <tumor_dna_id>.esvee.raw.vcf.gz.tbi            # Tumor DNA sample: VCF with candidate breakends
├── depth_annotation
│   ├── <tumor_dna_id>.esvee.ref_depth.vcf.gz          # Tumor DNA sample: VCF annotated with normal sample read depths
│   └── <tumor_dna_id>.esvee.ref_depth.vcf.gz.tbi      # Tumor DNA sample: VCF index
└── caller
    ├── <tumor_dna_id>.esvee.germline.vcf.gz           # Tumor DNA sample: VCF with germline breakends
    ├── <tumor_dna_id>.esvee.germline.vcf.gz.tbi       # Tumor DNA sample: VCF index
    ├── <tumor_dna_id>.esvee.somatic.vcf.gz            # Tumor DNA sample: VCF with somatic breakends
    ├── <tumor_dna_id>.esvee.somatic.vcf.gz.tbi        # Tumor DNA sample: VCF index
    ├── <tumor_dna_id>.esvee.unfiltered.vcf.gz         # Tumor DNA sample: VCF with unfiltered breakends
    └── <tumor_dna_id>.esvee.unfiltered.vcf.gz.tbi     # Tumor DNA sample: VCF index
```

### CNV calling

#### AMBER

[AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber) generates B-allele frequencies in tumor samples
for CNV calling in PURPLE.

```shell
<group_id>/amber/
├── <tumor_dna_id>.amber.baf.pcf                   # Tumor DNA sample: Piecewise constant fit on B-allele frequencies
├── <tumor_dna_id>.amber.baf.tsv.gz                # Tumor DNA sample: B-allele frequencies
├── <tumor_dna_id>.amber.contamination.tsv         # Tumor DNA sample: Contamination TSV
├── <tumor_dna_id>.amber.contamination.vcf.gz      # Tumor DNA sample: Contamination sites
├── <tumor_dna_id>.amber.contamination.vcf.gz.tbi  # Tumor DNA sample: Sample contamination sites index
├── <tumor_dna_id>.amber.qc                        # Tumor DNA sample: QC file
├── <normal_dna_id>.amber.homozygousregion.tsv     # Normal DNA sample: Regions of homozygosity
├── <normal_dna_id>.amber.snp.vcf.gz               # Normal DNA sample: SNP sites VCF
├── <normal_dna_id>.amber.snp.vcf.gz.tbi           # Normal DNA sample: VCF index
└── amber.version                                  # Tool version
```

#### COBALT

[COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt) generates read depth ratios (or an estimation
for tumor-only) for CNV calling in PURPLE.

```shell
<group_id>/cobalt/
├── <tumor_dna_id>.cobalt.gc.median.tsv      # Tumor DNA sample: GC median read depths
├── <tumor_dna_id>.cobalt.ratio.pcf          # Tumor DNA sample: Piecewise constant fit
├── <tumor_dna_id>.cobalt.ratio.tsv.gz       # Tumor DNA sample: Read counts and ratios (with reference or supposed diploid)
├── <normal_dna_id>.cobalt.gc.median.tsv     # Normal DNA sample: GC median read depths
├── <normal_dna_id>.cobalt.ratio.median.tsv  # Normal DNA sample: Chromosome median ratios
├── <normal_dna_id>.cobalt.ratio.pcf         # Normal DNA sample: Piecewise constant fit
└── cobalt.version                           # Tool version
```

#### PURPLE

[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) is a CNV caller that also infers tumor
purity/ploidy and annotates both small and structural variant calls with copy-number information.

```shell
<group_id>/purple/
├── <tumor_dna_id>.purple.cnv.gene.tsv                 # Tumor DNA sample: Somatic gene copy number
├── <tumor_dna_id>.purple.cnv.somatic.tsv              # Tumor DNA sample: Copy number variant segments
├── <tumor_dna_id>.purple.driver.catalog.germline.tsv  # Tumor DNA sample: Germline DNA sample driver events
├── <tumor_dna_id>.purple.driver.catalog.somatic.tsv   # Tumor DNA sample: Somatic DNA sample driver events
├── <tumor_dna_id>.purple.germline.deletion.tsv        # Tumor DNA sample: Germline DNA deletions
├── <tumor_dna_id>.purple.germline.vcf.gz              # Tumor DNA sample: Germline SAGE SMNVs with PURPLE annotations
├── <tumor_dna_id>.purple.germline.vcf.gz.tbi          # Tumor DNA sample: VCF index
├── <tumor_dna_id>.purple.purity.range.tsv             # Tumor DNA sample: Purity/ploidy model fit scores across a range of purity values
├── <tumor_dna_id>.purple.purity.tsv                   # Tumor DNA sample: Purity/ploidy summary
├── <tumor_dna_id>.purple.qc                           # Tumor DNA sample: QC file
├── <tumor_dna_id>.purple.segment.tsv                  # Tumor DNA sample: Genomic copy number segments
├── <tumor_dna_id>.purple.somatic.clonality.tsv        # Tumor DNA sample: Clonality peak model data
├── <tumor_dna_id>.purple.somatic.hist.tsv             # Tumor DNA sample: Somatic variants histogram data
├── <tumor_dna_id>.purple.somatic.vcf.gz               # Tumor DNA sample: Tumor SAGE SMNVs with PURPLE annotations
├── <tumor_dna_id>.purple.somatic.vcf.gz.tbi           # Tumor DNA sample: VCF index
├── <tumor_dna_id>.purple.sv.germline.vcf.gz           # Tumor DNA sample: Germline ESVEE SVs with PURPLE annotations
├── <tumor_dna_id>.purple.sv.germline.vcf.gz.tbi       # Tumor DNA sample: VCF index
├── <tumor_dna_id>.purple.sv.vcf.gz                    # Tumor DNA sample: Somatic ESVEE SVs with PURPLE annotations
├── <tumor_dna_id>.purple.sv.vcf.gz.tbi                # Tumor DNA sample: VCF index
├── circos/                                            # Circos plot data
├── plot/                                              # PURPLE plots
└── purple.version                                     # Tool version
```

### SV and driver event interpretation

#### LINX

[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx) clusters PURPLE-annotated SVs into high-order events
and classifies these events within a biological context. Following clustering and interpretation, events are visualised
as LINX plots.

```shell
<group_id>/linx/
├── germline_annotations
│   ├── <tumor_dna_id>.linx.germline.breakend.tsv        # Normal DNA sample: Breakend data
│   ├── <tumor_dna_id>.linx.germline.clusters.tsv        # Normal DNA sample: Clustered events
│   ├── <tumor_dna_id>.linx.germline.disruption.tsv      # Normal DNA sample: Structural disruptions
│   ├── <tumor_dna_id>.linx.germline.driver.catalog.tsv  # Normal DNA sample: Driver events
│   ├── <tumor_dna_id>.linx.germline.links.tsv           # Normal DNA sample: Cluster links
│   ├── <tumor_dna_id>.linx.germline.svs.tsv             # Normal DNA sample: Structural variants
│   └── linx.version                                     # Tool version
├── somatic_annotations
│   ├── <tumor_dna_id>.linx.breakend.tsv                 # Tumor DNA sample: Breakend data
│   ├── <tumor_dna_id>.linx.clusters.tsv                 # Tumor DNA sample: Clustered events
│   ├── <tumor_dna_id>.linx.driver.catalog.tsv           # Tumor DNA sample: Driver events
│   ├── <tumor_dna_id>.linx.drivers.tsv                  # Tumor DNA sample: Driver catalog
│   ├── <tumor_dna_id>.linx.fusion.tsv                   # Tumor DNA sample: Fusions
│   ├── <tumor_dna_id>.linx.links.tsv                    # Tumor DNA sample: Cluster links
│   ├── <tumor_dna_id>.linx.neoepitope.tsv               # Tumor DNA sample: Neoepitopes
│   ├── <tumor_dna_id>.linx.svs.tsv                      # Tumor DNA sample: Structural variants
│   ├── <tumor_dna_id>.linx.vis_copy_number.tsv          # Tumor DNA sample: Visualization: copy number
│   ├── <tumor_dna_id>.linx.vis_fusion.tsv               # Tumor DNA sample: Visualization: fusions
│   ├── <tumor_dna_id>.linx.vis_gene_exon.tsv            # Tumor DNA sample: Visualization: gene exons
│   ├── <tumor_dna_id>.linx.vis_protein_domain.tsv       # Tumor DNA sample: Visualization: protein domains
│   ├── <tumor_dna_id>.linx.vis_segments.tsv             # Tumor DNA sample: Visualization: segments
│   ├── <tumor_dna_id>.linx.vis_sv_data.tsv              # Tumor DNA sample: Visualization: structural variants
│   └── linx.version
└── somatic_plots
    ├── all
    │   └── <tumor_dna_id>.*.png                         # All cluster plots
    └── reportable
        └── <tumor_dna_id>.*.png                         # Driver cluster plots
```

### RNA transcript analysis

#### ISOFOX

[Isofox](https://github.com/hartwigmedical/hmftools/tree/master/isofox) analyses RNA alignment data to quantify
transcripts, identify novel splice junctions, and call fusions.

```shell
<group_id>/isofox/
├── <tumor_rna_id>.isf.alt_splice_junc.csv  # Tumor RNA sample: Alternative splice junctions
├── <tumor_rna_id>.isf.fusions.csv          # Tumor RNA sample: Fusions, unfiltered
├── <tumor_rna_id>.isf.gene_collection.csv  # Tumor RNA sample: Gene-collection fragment counts
├── <tumor_rna_id>.isf.gene_data.csv        # Tumor RNA sample: Gene fragment counts
├── <tumor_rna_id>.isf.pass_fusions.csv     # Tumor RNA sample: Fusions, filtered
├── <tumor_rna_id>.isf.retained_intron.csv  # Tumor RNA sample: Retained introns
├── <tumor_rna_id>.isf.summary.csv          # Tumor RNA sample: Analysis summary
└── <tumor_rna_id>.isf.transcript_data.csv  # Tumor RNA sample: Transcript fragment counts
```

### Oncoviral detection

#### VIRUSBreakend

[VIRUSBreakend](VIRUSBreakend) detects the presence of oncoviruses and intergration sites in tumor samples.

```shell
<group_id>/virusbreakend/
├── <tumor_dna_id>.virusbreakend.vcf              # Tumor DNA sample: VCF with viral integration sites
└── <tumor_dna_id>.virusbreakend.vcf.summary.tsv  # Tumor DNA sample: Analysis summary
```

#### VirusInterpreter

[VirusInterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter) performs post-processing for
VIRUSBreakend calls, such as merging of similar virus species into groups.

```shell
<group_id>/virusinterpreter/
└── <tumor_dna_id>.virus.annotated.tsv  # Tumor DNA sample: Processed oncoviral call/annotation data
```

### Immune analysis

#### LILAC

[LILAC](https://github.com/hartwigmedical/hmftools/tree/master/lilac) calls HLA Class I and characterises allelic status
(copy-number alterations, somatic mutations) in the tumor sample. Analysis can also incorporate RNA data as an
indirect measurement of allele expression.


```shell
<group_id>/lilac/
├── <tumor_dna_id>.lilac.candidates.coverage.tsv  # Tumor DNA sample: Coverage of high scoring candidates
├── <tumor_dna_id>.lilac.qc.tsv                   # Tumor DNA sample: QC file
└── <tumor_dna_id>.lilac.tsv                      # Tumor DNA sample: Analysis summary
```

#### NEO

[NEO](https://github.com/hartwigmedical/hmftools/tree/master/neo) identifies neoepitopes from point mutations, small indels and gene fusions,
as well as calculating allele specific neoepitope binding and presentation likelihoods.

```shell
<group_id>/neo/
├── <tumor_dna_id>.lilac.candidates.coverage.tsv  # Tumor DNA sample: Coverage of high scoring candidates
├── <tumor_dna_id>.lilac.qc.tsv                   # Tumor DNA sample: QC file
└── <tumor_dna_id>.lilac.tsv                      # Tumor DNA sample: Analysis summary
```

### Mutational signature fitting

#### SIGS

[Sigs](https://github.com/hartwigmedical/hmftools/tree/master/sigs) fits defined COSMIC trinucleotide mutational
signatures to tumor sample data.

```shell
sigs/
├── <tumor_dna_id>.sig.allocation.tsv  # Tumor DNA sample: SNV counts fitted to mutational signatures
└── <tumor_dna_id>.sig.snv_counts.csv  # Tumor DNA sample: SNV counts per 96 trinucleotide context
```

### HRD prediction

#### CHORD

[CHORD](https://github.com/UMCUGenetics/CHORD) is a random forest classifier of HRD status using relative counts of somatic mutations as
features.


```shell
<group_id>/chord/
├── <tumor_dna_id>.chord.mutation_contexts.tsv  # Tumor DNA sample: Counts of mutation types
└── <tumor_dna_id>.chord.prediction.tsv         # Tumor DNA sample: HRD predictions
```

### Tissue of origin prediction

#### CUPPA

[CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) predicts tissue of origin for a given tumor sample
using DNA and/or RNA features generated by upstream hmftools components.

```shell
<group_id>/cuppa/
├── <tumor_dna_id>.cuppa.pred_summ.tsv  # Tumor DNA sample: Prediction summary
├── <tumor_dna_id>.cuppa.vis.png        # Tumor DNA sample: Prediction visualisation
├── <tumor_dna_id>.cuppa.vis_data.tsv   # Tumor DNA sample: Prediction visualisation raw data
└── <tumor_dna_id>.cuppa_data.tsv.gz    # Tumor DNA sample: Input features
```

### Summary report

#### ORANGE

[ORANGE](https://github.com/hartwigmedical/hmftools/tree/master/orange) summaries and integrates key results from
`WiGiTS` components into a single static PDF report.

```shell
<group_id>/orange/
├── <tumor_dna_id>.orange.pdf   # Tumor DNA sample: Results of all tools as a PDF
└── <tumor_dna_id>.orange.json  # Tumor DNA sample: Result raw data
```

### Pipeline information

#### Created by Nextflow
```shell
pipeline_info/
├── execution_report_<date_time>.html    # HTML report of execution metrics and details
├── execution_timeline_<date_time>.html  # Timeline diagram showing process start/duration/finish
├── execution_trace_<date_time>.txt      # Resource usage
├── pipeline_dag_<date_time>.html        # Pipeline diagram showing how each process is connected
```

#### Created by Oncoanalyser
```shell
├── params_<date_time>.json              # Parameters used by the pipeline run
├── pipeline_report_<date_time>.html     # Only present if --email / --email_on_fail args are used
├── pipeline_report_<date_time>.txt      # Only present if --email / --email_on_fail args are used
└── software_versions.yml                # Tool versions
```
