# nf-core/oncoanalyser: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below will be created in the results
directory after the pipeline has finished. All paths are relative to the top-level results directory.

```text
output/
│  
├── subject_1/
│   ├── alignments/
│   ├── amber/
│   ├── bamtools/
│   ├── chord/
│   ├── cobalt/
│   ├── cuppa/
│   ├── flagstats/
│   ├── gridss/
│   ├── gripss/
│   ├── isofox/
│   ├── lilac/
│   ├── linx/
│   ├── orange/
│   ├── pave/
│   ├── purple/
│   ├── sage/
│   ├── sigs/
│   ├── virusbreakend/
│   └── virusinterpreter/
│  
├── subject_2/
│   └── ...
│  
...
│  
└── pipeline_info/
```

## Pipeline overview

- [Simple DNA/RNA alignment](#simple-dnarna-alignment)
  - [bwa-mem2](#bwa-mem2) - DNA alignment
  - [STAR](#star) - RNA alignment
- [Alignment post-processing](#alignment-post-processing)
  - [MarkDups](#markdups) - General alignment processing
  - [Picard Markduplicates](#picard-markduplicates) - Duplicate read marking
- [SNV, MNV, INDEL calling](#snv-mnv-indel-calling)
  - [SAGE](#sage) - SNV, MNV, INDEL calling
  - [PAVE](#pave) - Small variant annotation (transcript/coding effects)
- [SV calling](#sv-calling)
  - [SvPrep](#svprep) - Read filtering for SV calling
  - [GRIDSS](#gridss) - SV calling
  - [GRIPSS](#gripss) - SV filtering and post-processing
- [CNV calling](#cnv-calling)
  - [AMBER](#amber) - β-allele frequencies
  - [COBALT](#cobalt) - Read depth ratios
  - [PURPLE](#purple) - Purity/ploid estimation, variant annotation
- [SV event interpretation](#sv-event-interpretation)
  - [LINX](#linx) - SV event clustering and annotation
- [Transcript analysis](#transcript-analysis)
  - [Isofox](#isofox) - transcript counts, novel splicing and fusion calling
- [Oncoviral detection](#oncoviral-detection)
  - [VIRUSBreakend](#virusbreakend) - viral content and integration calling
  - [Virus Interpreter](#virus-interpreter) - oncoviral calling post-processing
- [HLA calling](#hla-calling)
  - [LILAC](#lilac) - HLA calling
- [HRD status prediction](#hrd-status-prediction)
  - [CHORD](#chord) - HRD status prediction
- [Mutational signature fitting](#mutational-signature-fitting)
  - [Sigs](#sigs) - Mutational signature fitting
- [Tissue of origin prediction](#tissue-of-origin-prediction)
  - [CUPPA](#cuppa) - Tissue of origin prediction
- [Report generation](#report-generation)
  - [ORANGE](#orange) - Key results summary
  - [linxreport](#linxreport) - Interactive LINX report
- [Pipeline information](#pipeline-information) - Workflow execution metrics

### Simple DNA/RNA alignment

Alignment functionality in oncoanalyser is simple and rigid, and exists only to meet the exact requirements of the
hmftools.

#### bwa-mem2

[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) is a short-read mapping tool used to align reads to a large reference
sequences. In oncoanalyser, bwa-mem2 is used to align DNA reads to the human genome.

_No outputs are published directly from bwa-mem2, see [MarkDups](#markdups) for the fully processed alignment outputs_

#### STAR

[STAR](https://github.com/alexdobin/STAR) is a specialised mapping to used to align RNA reads to a reference
transcriptome.

_No outputs are published directly from STAR, see [Picard MarkDuplicates](#picard-markduplicates) for the fully processed alignment outputs_

### Alignment post-processing

#### MarkDups

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/alignments/dna/`
  - `<normal_dna_id>.duplicate_freq.tsv`: Normal DNA sample read duplicate frequencies.
  - `<normal_dna_id>.markdups.bam`: Normal DNA sample output read alignments.
  - `<normal_dna_id>.markdups.bam.bai`: Normal DNA sample output read alignments index.
  - `<tumor_rna_id>.duplicate_freq.tsv`: Tumor DNA sample read duplicate frequencies.
  - `<tumor_rna_id>.markdups.bam`: Tumor DNA sample output read alignments.
  - `<tumor_rna_id>.markdups.bam.bai`: Tumor DNA sample output read alignments index.

</details>

[MarkDups](https://github.com/hartwigmedical/hmftools/tree/master/mark-dups) applies various alignment post-processing
routines such as duplicate marking and unmapping of problematic regions. It can also handle UMIs when configured to do
so.

_MarkDups is only run on DNA alignments_

### Picard MarkDuplicates

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/alignments/rna/`
  - `<tumor_rna_id>.md.bam`: Tumor RNA sample read alignments.
  - `<tumor_rna_id>.md.bam.bai`: Tumor RNA sample read alignments index.
  - `<tumor_rna_id>.md.metrics`: Tumor RNA sample read duplicate marking metrics.

</details>

[Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) used to
mark duplicate reads following alignment.

_Picard MarkDuplicates is only run on RNA alignments_

### SNV, MNV, INDEL calling

#### SAGE

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/sage/append/`

  - `<tumor_dna_id>.sage.append.vcf.gz`: Tumor DNA sample small variant VCF with RNA data appended.
  - `<normal_dna_id>.sage.append.vcf.gz`: Normal DNA sample small variant VCF with RNA data appended.

- `<group_id>/sage/somatic/`

  - `<normal_dna_id>.sage.bqr.png`: Normal DNA sample base quality recalibration metrics plot.
  - `<normal_dna_id>.sage.bqr.tsv`: Normal DNA sample base quality recalibration metrics.
  - `<tumor_dna_id>.sage.bqr.png`: Tumor DNA sample base quality recalibration metrics plot.
  - `<tumor_dna_id>.sage.bqr.tsv`: Tumor DNA sample base quality recalibration metrics.
  - `<tumor_dna_id>.sage.exon.medians.tsv`: Tumor DNA sample exon median depths.
  - `<tumor_dna_id>.sage.gene.coverage.tsv`: Tumor DNA sample gene coverages.
  - `<tumor_dna_id>.sage.somatic.filtered.vcf.gz.tbi`: Tumor DNA sample filtered small variant calls index.
  - `<tumor_dna_id>.sage.somatic.filtered.vcf.gz`: Tumor DNA sample filtered small variant calls.
  - `<tumor_dna_id>.sage.somatic.vcf.gz.tbi`: Tumor DNA sample small variant calls index.
  - `<tumor_dna_id>.sage.somatic.vcf.gz`: Tumor DNA sample small variant calls.

- `<group_id>/sage/germline/`
  - `<normal_dna_id>.sage.bqr.png`: Tumor DNA sample base quality recalibration metrics plot.
  - `<normal_dna_id>.sage.bqr.tsv`: Tumor DNA sample base quality recalibration metrics.
  - `<normal_dna_id>.sage.exon.medians.tsv`: Normal DNA sample exon median depths.
  - `<normal_dna_id>.sage.gene.coverage.tsv`: Normal DNA sample gene coverages.
  - `<tumor_dna_id>.sage.bqr.png`: Normal DNA sample base quality recalibration metrics plot.
  - `<tumor_dna_id>.sage.bqr.tsv`: Normal DNA sample base quality recalibration metrics.
  - `<tumor_dna_id>.sage.germline.filtered.vcf.gz.tbi`: Normal DNA sample filtered small variant calls index.
  - `<tumor_dna_id>.sage.germline.filtered.vcf.gz`: Normal DNA sample filtered small variant calls.
  - `<tumor_dna_id>.sage.germline.vcf.gz.tbi`: Normal DNA sample small variant calls index.
  - `<tumor_dna_id>.sage.germline.vcf.gz`: Normal DNA sample small variant calls.

</details>

[SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage) is a SNV, MNV, and INDEL caller optimised for 100x
tumor and 40x normal.

#### PAVE

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/pave/`
  - `<tumor_dna_id>.sage.germline.filtered.pave.vcf.gz.tbi`: Annotated SAGE germline small variants index.
  - `<tumor_dna_id>.sage.germline.filtered.pave.vcf.gz`: Annotated SAGE germline small variants.
  - `<tumor_dna_id>.sage.somatic.filtered.pave.vcf.gz.tbi`: Annotated SAGE somatic small variants index.
  - `<tumor_dna_id>.sage.somatic.filtered.pave.vcf.gz`: Annotated SAGE somatic small variants.

</details>

[PAVE](https://github.com/hartwigmedical/hmftools/tree/master/pave) annotates variants called by SAGE with impact
information with regards to transcript and coding effects.

### SV calling

#### SvPrep

[SvPrep](https://github.com/hartwigmedical/hmftools/tree/master/sv-prep) runs prior to SV calling to reducing runtime by
rapidly identifying reads that are likely to be involved in a SV event.

_No outputs are published directly from SvPrep, see [GRIPSS](#gripss) for the fully processed SV calling outputs_

#### GRIDSS

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/gridss/`
  - `<tumor_dna_id>.gridss.vcf.gz`: GRIDSS structural variants.
  - `<tumor_dna_id>.gridss.vcf.gz.tbi`: GRIDSS structural variants index.

</details>

[GRIDSS](https://github.com/PapenfussLab/gridss) is a SV caller than uses both read support and local
breakend/breakpoint assemblies to call variants.

#### GRIPSS

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/gripss/germline/`

  - `<tumor_dna_id>.gripss.filtered.germline.vcf.gz`: Filtered GRIDSS germline structural variants.
  - `<tumor_dna_id>.gripss.filtered.germline.vcf.gz.tbi`: Filtered GRIDSS germline structural variants index.
  - `<tumor_dna_id>.gripss.germline.vcf.gz`: GRIDSS structural variants (GRIPSS filters set but not applied).
  - `<tumor_dna_id>.gripss.germline.vcf.gz.tbi`: GRIDSS structural variants index (GRIPSS filters set but not applied).

- `<group_id>/gripss/somatic/`
  - `<tumor_dna_id>.gripss.filtered.somatic.vcf.gz`: Filtered GRIDSS somatic structural variants.
  - `<tumor_dna_id>.gripss.filtered.somatic.vcf.gz.tbi`: Filtered GRIDSS somatic structural variants index.
  - `<tumor_dna_id>.gripss.somatic.vcf.gz`: GRIDSS structural variants (GRIPSS filters set but not applied).
  - `<tumor_dna_id>.gripss.somatic.vcf.gz.tbi`: GRIDSS structural variants index (GRIPSS filters set but not applied).

</details>

[GRIPSS](https://github.com/hartwigmedical/hmftools/tree/master/gripss) applies filter and post-processing to SV calls.

### CNV calling

#### AMBER

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/amber/`
  - `amber.version`: AMBER version file.
  - `<tumor_dna_id>.amber.baf.pcf`: Tumor DNA sample piecewise constant fit.
  - `<tumor_dna_id>.amber.baf.tsv.gz`: Tumor DNA sample β-allele frequencies.
  - `<tumor_dna_id>.amber.contamination.tsv`: Tumor DNA sample contamination TSV.
  - `<tumor_dna_id>.amber.contamination.vcf.gz`: Tumor DNA sample contamination sites.
  - `<tumor_dna_id>.amber.contamination.vcf.gz.tbi`: Tumor DNA sample contamination sites index.
  - `<tumor_dna_id>.amber.qc`: AMBER QC file.
  - `<normal_dna_id>.amber.homozygousregion.tsv`: Normal DNA sample regions of homozygosity.

</details>

[AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber) generates β-allele frequencies in tumor samples
for CNV calling in PURPLE.

#### COBALT

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/cobalt/`
  - `cobalt.version`: COBALT version file.
  - `<tumor_dna_id>.cobalt.gc.median.tsv`: Tumor DNA sample GC median read depths.
  - `<tumor_dna_id>.cobalt.ratio.pcf`: Tumor DNA sample piecewise constant fit.
  - `<tumor_dna_id>.cobalt.ratio.tsv.gz`: Tumor DNA sample read counts and ratios (with reference or supposed diploid
    regions).
  - `<normal_dna_id>.cobalt.gc.median.tsv`: Normal DNA sample GC median read depths.
  - `<normal_dna_id>.cobalt.ratio.median.tsv`: Normal DNA sample chromosome median ratios.
  - `<normal_dna_id>.cobalt.ratio.pcf`: Normal DNA sample piecewise constant fit.

</details>

[COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt) generates read depth ratios (or an estimation
for tumor-only) for CNV calling in PURPLE.

#### PURPLE

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/purple/`
  - `circos/`: Circos plot data.
  - `<tumor_dna_id>.purple.cnv.gene.tsv`: Somatic gene copy number.
  - `<tumor_dna_id>.purple.cnv.somatic.tsv`: Copy number variant segments.
  - `<tumor_dna_id>.purple.driver.catalog.germline.tsv`: Normal DNA sample driver catalogue.
  - `<tumor_dna_id>.purple.driver.catalog.somatic.tsv`: Tumor DNA sample driver catalogue.
  - `<tumor_dna_id>.purple.germline.deletion.tsv`: Normal DNA deletions.
  - `<tumor_dna_id>.purple.germline.vcf.gz`: Normal DNA SAGE small variants with PURPLE annotations.
  - `<tumor_dna_id>.purple.germline.vcf.gz.tbi`: Normal DNA SAGE small variants with PURPLE annotations index.
  - `<tumor_dna_id>.purple.purity.range.tsv`: Purity/ploid model fit scores across a range of purity values.
  - `<tumor_dna_id>.purple.purity.tsv`: Purity/ploidy summary.
  - `<tumor_dna_id>.purple.qc`: PURPLE QC file.
  - `<tumor_dna_id>.purple.segment.tsv`: Genomic copy number segments.
  - `<tumor_dna_id>.purple.somatic.clonality.tsv`: Clonality peak model data.
  - `<tumor_dna_id>.purple.somatic.hist.tsv`: Somatic variants histogram data.
  - `<tumor_dna_id>.purple.somatic.vcf.gz`: Tumor DNA sample small variants with PURPLE annotations.
  - `<tumor_dna_id>.purple.somatic.vcf.gz.tbi`: Tumor DNA sample small variants with PURPLE annotations index.
  - `<tumor_dna_id>.purple.sv.germline.vcf.gz`: Germline structural variants with PURPLE annotations.
  - `<tumor_dna_id>.purple.sv.germline.vcf.gz.tbi`: Germline structural variants with PURPLE annotations index.
  - `<tumor_dna_id>.purple.sv.vcf.gz`: Somatic structural variants with PURPLE annotations.
  - `<tumor_dna_id>.purple.sv.vcf.gz.tbi`: Somatic structural variants with PURPLE annotations.
  - `plot/`: PURPLE plots.
  - `purple.version`: PURPLE version file.

</details>

[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) is a CNV caller that also infers tumor
purity/ploidy and annotates both small and structural variant calls with copy-number information.

### SV event interpretation

#### LINX

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/linx/germline_annotations/`

  - `linx.version`: LINX version file.
  - `<tumor_dna_id>.linx.germline.breakend.tsv`: Normal DNA sample breakend data.
  - `<tumor_dna_id>.linx.germline.clusters.tsv`: Normal DNA sample clustered events.
  - `<tumor_dna_id>.linx.germline.disruption.tsv`: Normal DNA sample breakend data.
  - `<tumor_dna_id>.linx.germline.driver.catalog.tsv`: Normal DNA sample driver catalogue.
  - `<tumor_dna_id>.linx.germline.links.tsv`: Normal DNA sample cluster links.
  - `<tumor_dna_id>.linx.germline.svs.tsv`: Normal DNA sample structural variants.

- `<group_id>/linx/somatic_annotations/`

  - `linx.version`: LINX version file.
  - `<tumor_dna_id>.linx.breakend.tsv`: Tumor DNA sample breakend data.
  - `<tumor_dna_id>.linx.clusters.tsv`: Tumor DNA sample clustered events.
  - `<tumor_dna_id>.linx.driver.catalog.tsv`: Tumor DNA sample driver catalogue.
  - `<tumor_dna_id>.linx.drivers.tsv`: Tumor DNA sample LINX driver drivers.
  - `<tumor_dna_id>.linx.fusion.tsv`: Tumor DNA sample fusions.
  - `<tumor_dna_id>.linx.links.tsv`: Tumor DNA sample cluster links.
  - `<tumor_dna_id>.linx.svs.tsv`: Tumor DNA sample structural variants.
  - `<tumor_dna_id>.linx.vis_*`: Tumor DNA sample visualisation data.

- `<group_id>/linx/somatic_plots/`
  - `all/*png`: All available tumor DNA sample cluster plots.
  - `reportable/*png`: Driver-only tumor DNA sample cluster plots.

</details>

[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx) clusters PURPLE-annotated SVs into high-order events
and classifies these events within a biological context. Following clustering and interpretation, events are visualised
as LINX plots.

### Transcript analysis

#### Isofox

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/isofox/`
  - `<tumor_rna_id>.isf.alt_splice_junc.csv`: Tumor RNA sample alternative splice junctions.
  - `<tumor_rna_id>.isf.fusions.csv`: Tumor RNA sample fusions, unfiltered.
  - `<tumor_rna_id>.isf.gene_collection.csv`: Tumor RNA sample gene-collection fragment counts.
  - `<tumor_rna_id>.isf.gene_data.csv`: Tumor RNA sample gene fragment counts.
  - `<tumor_rna_id>.isf.pass_fusions.csv`: Tumor RNA sample fusions, filtered.
  - `<tumor_rna_id>.isf.retained_intron.csv`: Tumor RNA sample retained introns.
  - `<tumor_rna_id>.isf.summary.csv`: Tumor RNA sample analysis summary file.
  - `<tumor_rna_id>.isf.transcript_data.csv`: Tumor RNA sample transcript fragment counts.

</details>

[Isofox](https://github.com/hartwigmedical/hmftools/tree/master/isofox) analyses RNA alignment data to quantify
transcripts, identify novel splice junctions, and caller fusions.

### Oncoviral detection

#### VIRUSBreakend

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/virusbreakend/`
  - `<tumor_dna_id>.virusbreakend.vcf`: Tumor DNA sample viral integratino sites.
  - `<tumor_dna_id>.virusbreakend.vcf.summary.tsv`: Tumor DNA sample analysis summary file.

</details>

[VIRUSBreakend](VIRUSBreakend) detects the presence of oncoviruses and intergration sites in tumor samples.

#### Virus Interpreter

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/virusinterpreter/`
  - `<tumor_dna_id>.virus.annotated.tsv`: Processed oncoviral call/annotation data.

</details>

[Virus Interpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter) post-processing for
VIRUSBreakend calls that provides higher-level interpretation of data.

### HLA calling

#### LILAC

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/lilac/`
  - `<tumor_dna_id>.lilac.candidates.coverage.tsv`: Coverage of high scoring candidates.
  - `<tumor_dna_id>.lilac.qc.tsv`: LILAC qc file.
  - `<tumor_dna_id>.lilac.tsv`: Analysis summary.

</details>

[LILAC](https://github.com/hartwigmedical/hmftools/tree/master/lilac) calls HLA Class I and characterises allelic status
(copy-number alterations, somatic mutations) in the tumor sample. Analysis can also incorporate RNA data as an
indirectly measurement of allele expression.

### HRD status prediction

#### CHORD

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/chord/`
  - `<tumor_dna_id>_chord_prediction.txt`: Tumor DNA sample analysis summary file.
  - `<tumor_dna_id>_chord_signatures.txt`: Tumor DNA sample variant counts contributing to signatures.

</details>

[CHORD](https://github.com/UMCUGenetics/CHORD) predicts the HRD status of a tumor using statistical inference on the
basis of relative somatic mutation counts.

### Mutational signature fitting

#### Sigs

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/sigs/`
  - `<tumor_dna_id>.sig.allocation.tsv`: Tumor DNA sample signature allocations.
  - `<tumor_dna_id>.sig.snv_counts.csv`: Tumor DNA sample variant counts contributing to signatures.

</details>

[Sigs](https://github.com/hartwigmedical/hmftools/tree/master/sigs) fits defined COSMIC trinucleotide mutational
signatures to tumor sample data.

### Tissue of origin prediction

#### CUPPA

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/cuppa/`
  - `<tumor_dna_id>_cup_report.pdf`: Combined figure of summary and feature plot.
  - `<tumor_dna_id>.cup.data.csv`: Model feature scores.
  - `<tumor_dna_id>.cup.report.features.png`: Feature plot.
  - `<tumor_dna_id>.cup.report.summary.png`: Summary plot.
  - `<tumor_dna_id>.cuppa.chart.png`: CUPPA chart plot.
  - `<tumor_dna_id>.cuppa.conclusion.txt`: Prediction conclusion file.

</details>

[CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) predicts tissue of origin for a given tumor sample
using DNA and/or RNA features generated by upstream hmftools components.

### Report generation

#### ORANGE

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/orange/`
  - `<tumor_dna_id>.orange.json`: Aggregated report data.
  - `<tumor_dna_id>.orange.pdf`: Static report PDF.

</details>

[ORANGE](https://github.com/hartwigmedical/hmftools/tree/master/orange) summaries and integrates key results from
hmftool components into a single static PDF report.

#### linxreport

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/linx/`
  - `<tumor_dna_id>_linx.html`: Interactive HTML report.

</details>

[linxreport](https://github.com/umccr/linxreport) generates an interactive report containing LINX annotations and plots.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>
