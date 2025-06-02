# nf-core/oncoanalyser: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below will be created in the results
directory after the pipeline has finished. All paths are relative to the top-level results directory.

```tree
output/
│  
├── <group_id-1>/
│   ├── alignments/
│   ├── amber/
│   ├── bamtools/
│   ├── chord/
│   ├── cider/
│   ├── cobalt/
│   ├── cuppa/
│   ├── esvee/
│   ├── isofox/
│   ├── lilac/
│   ├── linx/
│   ├── orange/
│   ├── pave/
│   ├── peach/
│   ├── purple/
│   ├── sage/
│   ├── sigs/
│   ├── teal/
│   ├── virusbreakend/
│   └── virusinterpreter/
│  
├── <group_id-2>/
│   └── ...
│  
...
│  
└── pipeline_info/
```

## Pipeline overview

- [Read alignment](#read-alignment)
  - [BWA-MEM2](#bwa-mem2) - DNA read alignment
  - [STAR](#star) - RNA read alignment
- [Read post-processing](#alignment-post-processing)
  - [REDUX](#redux) - DNA read post-processing
  - [Picard MarkDuplicates](#picard-markduplicates) - RNA read post-processing
- [SNV, MNV, INDEL calling](#snv-mnv-indel-calling)
  - [SAGE](#sage) - SNV, MNV, INDEL calling
  - [PAVE](#pave) - Small variant annotation (transcript/coding effects)
- [SV calling](#sv-calling)
  - [ESVEE](#esvee) - Read selection, SV calling, and variant filtering
- [CNV calling](#cnv-calling)
  - [AMBER](#amber) - β-allele frequencies
  - [COBALT](#cobalt) - Read depth ratios
  - [PURPLE](#purple) - Purity/ploid estimation, variant annotation
- [SV event interpretation](#sv-event-interpretation)
  - [LINX](#linx) - SV event clustering and annotation
- [Transcript analysis](#transcript-analysis)
  - [ISOFOX](#isofox) - RNA transcript analysis
- [Oncoviral detection](#oncoviral-detection)
  - [VIRUSBreakend](#virusbreakend) - Viral content and integration calling
  - [VirusInterpreter](#virusinterpreter) - Oncoviral calling post-processing
- [Telomere characterisation](#telomere-characterisation)
  - [TEAL](#teal) - Telomere characterisation
- [Immune analysis](#immune-analysis)
  - [LILAC](#lilac) - HLA typing
  - [CIDER](#cider) - IG/TCR CDR3 identification
  - [NEO](#neo) - Neoepitope prediction
- [Mutational signature fitting](#mutational-signature-fitting)
  - [SIGS](#sigs) - Mutational signature fitting
- [HRD status prediction](#hrd-status-prediction)
  - [CHORD](#chord) - HRD status prediction
- [Tissue of origin prediction](#tissue-of-origin-prediction)
  - [CUPPA](#cuppa) - Tissue of origin prediction
- [Pharmacogenomics](#pharmacogenomics)
  - [PEACH](#peach) - Pharmacogenomic assessment
- [Report generation](#report-generation)
  - [ORANGE](#orange) - Summary report
  - [linxreport](#linxreport) - Interactive LINX report
- [Pipeline information](#pipeline-information) - Workflow execution metrics

### Read alignment

Alignment functionality in `oncoanalyser` is simple and rigid, and exists only to meet the exact requirements of the
WiGiTS workflow.

#### BWA-MEM2

[BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) is a short-read mapping tool used to align reads to a large reference
sequences. In `oncoanalyser`, BWA-MEM2 is used to align DNA reads to the human genome.

_No outputs are published directly from bwa-mem2, see [REDUX](#redux) for the fully processed alignment outputs_

#### STAR

[STAR](https://github.com/alexdobin/STAR) is a specialised mapping tool used to align RNA reads to a reference
transcriptome.

_No outputs are published directly from STAR, see [Picard MarkDuplicates](#picard-markduplicates) for the fully processed alignment outputs_

### Alignment post-processing

#### REDUX

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/alignments/dna/`
  - `<tumor_dna_id>.jitter_params.tsv`: Tumor DNA sample microsatellite jitter model parameters.
  - `<tumor_dna_id>.ms_table.tsv.gz`: Tumor DNA sample aggregated repeat units and repeat counts.
  - `<tumor_dna_id>.redux.bam`: Tumor DNA sample read alignments.
  - `<tumor_dna_id>.redux.bam.bai`: Tumor DNA sample read alignments index.
  - `<tumor_dna_id>.redux.duplicate_freq.tsv`: Tumor DNA sample duplicate read frequencies.
  - `<tumor_dna_id>.repeat.tsv.gz`: Tumor DNA sample repeat units and repeat counts per site.
  - `<normal_dna_id>.jitter_params.tsv`: Normal DNA sample microsatellite jitter model parameters.
  - `<normal_dna_id>.ms_table.tsv.gz`: Normal DNA sample aggregated repeat units and repeat counts.
  - `<normal_dna_id>.redux.bam`: Normal DNA sample: Read alignments.
  - `<normal_dna_id>.redux.bam.bai`: Normal DNA sample read alignments index.
  - `<normal_dna_id>.redux.duplicate_freq.tsv`: Normal DNA sample duplicate read frequencies.
  - `<normal_dna_id>.repeat.tsv.gz`: Normal DNA sample repeat units and repeat counts per site.

</details>

[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) applies various alignment post-processing routines
such as duplicate marking and unmapping of problematic regions for DNA reads. It can also handle UMIs when configured to
do so.

_REDUX is only run on DNA alignments_

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
  - `<tumor_dna_id>.sage.somatic.vcf.gz.tbi`: Tumor DNA sample small variant calls index.
  - `<tumor_dna_id>.sage.somatic.vcf.gz`: Tumor DNA sample small variant calls.

- `<group_id>/sage/germline/`
  - `<normal_dna_id>.sage.bqr.png`: Normal DNA sample base quality recalibration metrics plot.
  - `<normal_dna_id>.sage.bqr.tsv`: Normal DNA sample base quality recalibration metrics.
  - `<normal_dna_id>.sage.exon.medians.tsv`: Normal DNA sample exon median depths.
  - `<normal_dna_id>.sage.gene.coverage.tsv`: Normal DNA sample gene coverages.
  - `<tumor_dna_id>.sage.bqr.png`: Normal DNA sample base quality recalibration metrics plot.
  - `<tumor_dna_id>.sage.bqr.tsv`: Normal DNA sample base quality recalibration metrics.
  - `<tumor_dna_id>.sage.germline.vcf.gz.tbi`: Normal DNA sample small variant calls index.
  - `<tumor_dna_id>.sage.germline.vcf.gz`: Normal DNA sample small variant calls.

</details>

[SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage) is a precise and highly sensitive somatic SNV, MNV
and small INDEL caller. It has dynamically scaling sensitivity based on the depth of the provided tumor and germline
BAMs, but performs best if both BAMs have at least 30x typical depth.

#### PAVE

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/pave/`
  - `<tumor_dna_id>.sage.germline.pave.vcf.gz.tbi`: Annotated SAGE germline small variants index.
  - `<tumor_dna_id>.sage.germline.pave.vcf.gz`: Annotated SAGE germline small variants.
  - `<tumor_dna_id>.sage.somatic.pave.vcf.gz.tbi`: Annotated SAGE somatic small variants index.
  - `<tumor_dna_id>.sage.somatic.pave.vcf.gz`: Annotated SAGE somatic small variants.

</details>

[PAVE](https://github.com/hartwigmedical/hmftools/tree/master/pave) annotates variants called by SAGE with impact
information with regards to transcript and coding effects.

### SV calling

#### ESVEE

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/esvee/prep/`

  - `<tumor_dna_id>.esvee.prep.bam`: Tumor DNA sample BAM with candidate SV reads.
  - `<tumor_dna_id>.esvee.prep.bam.bai`: Tumor DNA sample BAM index.
  - `<tumor_dna_id>.esvee.prep.disc_stats.tsv`: Tumor DNA sample discordant reads stats.
  - `<tumor_dna_id>.esvee.prep.fragment_length.tsv`: Tumor DNA sample fragment length stats.
  - `<tumor_dna_id>.esvee.prep.junction.tsv`: Tumor DNA sample candidate junctions.
  - `<normal_dna_id>.esvee.prep.bam`: Tumor DNA sample BAM with candidate SV reads.
  - `<normal_dna_id>.esvee.prep.bam.bai`: Tumor DNA sample BAM index.

- `<group_id>/esvee/assemble/`

  - `<tumor_dna_id>.esvee.assembly.tsv`: Tumor DNA sample breakend assemblies.
  - `<tumor_dna_id>.esvee.alignment.tsv`: Tumor DNA sample assemblies realigned to the reference genome.
  - `<tumor_dna_id>.esvee.breakend.tsv`: Tumor DNA sample breakends.
  - `<tumor_dna_id>.esvee.phased_assembly.tsv`: Tumor DNA sample phased assemblies.
  - `<tumor_dna_id>.esvee.raw.vcf.gz`: Tumor DNA sample VCF with candidate breakends.
  - `<tumor_dna_id>.esvee.raw.vcf.gz.tbi`: Tumor DNA sample VCF with candidate breakends.

- `<group_id>/esvee/depth_annotation/`

  - `<tumor_dna_id>.esvee.ref_depth.vcf.gz`: Tumor DNA sample VCF annotated with normal sample read depths.
  - `<tumor_dna_id>.esvee.ref_depth.vcf.gz.tbi`: Tumor DNA sample VCF index.

- `<group_id>/esvee/caller/`
  - `<tumor_dna_id>.esvee.germline.vcf.gz`: Tumor DNA sample VCF with germline breakends.
  - `<tumor_dna_id>.esvee.germline.vcf.gz.tbi`: Tumor DNA sample VCF index.
  - `<tumor_dna_id>.esvee.somatic.vcf.gz`: Tumor DNA sample VCF with somatic breakends.
  - `<tumor_dna_id>.esvee.somatic.vcf.gz.tbi`: Tumor DNA sample VCF index.
  - `<tumor_dna_id>.esvee.unfiltered.vcf.gz`: Tumor DNA sample VCF with unfiltered breakends.
  - `<tumor_dna_id>.esvee.unfiltered.vcf.gz.tbi`: Tumor DNA sample VCF index.

</details>

[ESVEE](https://github.com/hartwigmedical/hmftools/tree/master/esvee) is a SV caller that uses both read support and
local breakend/breakpoint assemblies to call variants.

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
  - `<normal_dna_id>.amber.snp.vcf.gz`: Normal DNA sample SNP sites.
  - `<normal_dna_id>.amber.snp.vcf.gz.tbi`: Normal DNA sample SNP sites index.

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
  - `<tumor_dna_id>.cobalt.ratio.tsv.gz`: Tumor DNA sample read counts and ratios (with reference or supposed diploid regions).
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
  - `<tumor_dna_id>.linx.drivers.tsv`: Tumor DNA sample LINX drivers.
  - `<tumor_dna_id>.linx.fusion.tsv`: Tumor DNA sample fusions.
  - `<tumor_dna_id>.linx.links.tsv`: Tumor DNA sample cluster links.
  - `<tumor_dna_id>.linx.neoepitope.tsv`: Tumor DNA sample neoepitopes.
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

#### ISOFOX

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

[ISOFOX](https://github.com/hartwigmedical/hmftools/tree/master/isofox) analyses RNA alignment data to quantify
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

#### VirusInterpreter

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/virusinterpreter/`
  - `<tumor_dna_id>.virus.annotated.tsv`: Processed oncoviral call/annotation data.

</details>

[VirusInterpreter](https://github.com/hartwigmedical/hmftools/tree/master/virus-interpreter) performs post-processing
for VIRUSBreakend calls to provide higher-level interpretation of data.

### Telomere characterisation

#### TEAL

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/teal/`
  - `<normal_dna_id>.teal.telbam.bam`: Normal DNA sample BAM with telomeric fragments.
  - `<normal_dna_id>.teal.telbam.bam.bai`: Normal DNA sample BAM index.
  - `<normal_dna_id>.teal.tellength.tsv`: Normal DNA sample telomeric content and length estimate.
  - `<tumor_dna_id>.teal.breakend.tsv.gz`: Tumor DNA sample telomeric rearrangements.
  - `<tumor_dna_id>.teal.telbam.bam`: Tumor DNA sample BAM with telomeric fragments.
  - `<tumor_dna_id>.teal.telbam.bam.bai`: Tumor DNA sample BAM index.
  - `<tumor_dna_id>.teal.tellength.tsv`: Tumor DNA sample telomeric content and length estimate.

</details>

[TEAL](https://github.com/hartwigmedical/hmftools/tree/master/teal) measures telomere content, and estimates telomeric
length based on WGS BAM input.

### Immune analysis

#### LILAC

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/lilac/`
  - `<tumor_dna_id>.lilac.candidates.coverage.tsv`: Coverage of high scoring candidates.
  - `<tumor_dna_id>.lilac.qc.tsv`: LILAC QC file.
  - `<tumor_dna_id>.lilac.tsv`: Analysis summary.

</details>

[LILAC](https://github.com/hartwigmedical/hmftools/tree/master/lilac) calls HLA Class I and characterises allelic status
(copy-number alterations, somatic mutations) in the tumor sample. Analysis can also incorporate RNA data as an indirect
measurement of allele expression.

#### CIDER

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/cider/`
  - `<tumor_dna_id|tumor_rna_id>.cider.bam`: Tumor DNA or RNA read alignments overlapping anchor sites.
  - `<tumor_dna_id|tumor_rna_id>.cider.blastn_match.tsv.gz`: Tumor DNA or RNA BLASTn hits of sequences.
  - `<tumor_dna_id|tumor_rna_id>.cider.layout.gz`: Tumor DNA or RNA read layout of each locus.
  - `<tumor_dna_id|tumor_rna_id>.cider.locus_stats.tsv`: Tumor DNA or RNA locus statistics and summary file.
  - `<tumor_dna_id|tumor_rna_id>.cider.vdj.tsv.gz`: Tumor DNA or RNA annotated VDJ sequences.

</details>

[CIDER](https://github.com/hartwigmedical/hmftools/tree/master/cider) determines a comprehensive list of CDR3 sequences
for each of the IG and TCR loci including an abundance estimate.

#### NEO

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/neo/finder/`

  - `<tumor_dna_id>.neo_data.tsv`: Tumor sample neoepitope candidates.

- `<group_id>/neo/annotated_fusions/`

  - `<tumor_dna_id>.isf.neoepitope.tsv`: Tumor sample annotated Isofox fusions.

- `<group_id>/neo/scorer/`
  - `<tumor_dna_id>.neo.peptide_scores.tsv`: Tumor sample peptide scores.
  - `<tumor_dna_id>.neo.neoepitope.tsv`: Tumor sample neoepitope predictions.

</details>

[NEO](https://github.com/hartwigmedical/hmftools/tree/master/neo) identifies neoepitopes from point mutations, small
indels and gene fusions, as well as calculating allele specific neoepitope binding and presentation likelihoods.

### HRD status prediction

#### CHORD

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/chord/`
  - `<tumor_dna_id>.chord.prediction.tsv`: Tumor DNA sample analysis summary file.
  - `<tumor_dna_id>.chord.mutation_contexts.tsv`: Tumor DNA sample counts of mutation types.

</details>

[CHORD](https://github.com/hartwigmedical/hmftools/tree/master/chord) predicts the HRD status of a tumor using
statistical inference on the basis of relative somatic mutation counts.

### Mutational signature fitting

#### SIGS

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
  - `<tumor_dna_id>.cuppa_data.tsv.gz`: Tumor sample input features
  - `<tumor_dna_id>.cuppa.pred_summ.tsv`: Tumor sample prediction summary
  - `<tumor_dna_id>.cuppa.vis_data.tsv`: Tumor sample predication visualisation data
  - `<tumor_dna_id>.cuppa.vis.png`: Tumor sample prediction visualisation

</details>

[CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) predicts tissue of origin for a given tumor sample
using DNA and/or RNA features generated by upstream WiGiTS components.

### Pharmacogenomics

#### PEACH

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/peach/`
  - `<normal_dna_id>.peach.events.tsv`: Normal DNA sample variant events.
  - `<normal_dna_id>.peach.gene.events.tsv`: Normal DNA sample variant events (linked by gene).
  - `<normal_dna_id>.peach.haplotypes.all.tsv`: Normal DNA sample all haplotypes.
  - `<normal_dna_id>.peach.haplotypes.best.tsv`: Normal DNA sample best haplotypes..
  - `<normal_dna_id>.peach.qc.tsv`: PEACH QC file.

</details>

[PEACH](https://github.com/hartwigmedical/hmftools/tree/master/peach) infers haplotypes for interpretation in a
pharmacogenomic context.

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

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
