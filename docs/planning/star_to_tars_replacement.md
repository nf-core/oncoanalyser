# TARS (Transcript Alignment for RNA Splicing) 

## What is TARS?

TARS makes bwa-mem2 splice-aware for RNA reads: it aligns against the genome with the transcriptome added, then rewrites the result back to genome coordinates. The output is an ordinary genomic RNA BAM (no transcript contigs, spliced reads carried as N gaps) ready for REDUX and ISOFOX. 

```mermaid
FASTQ_files --> BWA_MEM2 --> TARS --> REDUX --> ISOFOX
```

## Overview of what these 3 steps look like with ISOFOX at the end:

### 1. bwa-mem2

Align against the combined genome + transcript FASTA. Output must stay name-grouped, not coordinate sorted, because TARS consumes read groups

```bash
bwa-mem2 mem \
-R "@RG\tID:<lane>\tLB:<sample>\tPL:illumina\tPU:<sample>\tDS:Short\tSM:<sample>" \
-Y \
-T 19 \
-h 75 \
-t <threads> \
/path/to/Homo_sapiens_assembly38.alt.masked.with_rna_contigs.fasta \
<sample>_<lane>_R1_001.fastq.gz \
<sample>_<lane>_R2_001.fastq.gz \
| sambamba view -t <threads> -f bam -S /dev/stdin -o <sample>.<lane>.namegrouped.bam
```

Versions we have been running: bwa-mem2 20250204.7aa5ff6, sambamba 0.6.8

The sambamba step is only SAM to BAM, no sort, so happy to pick whatever we already have wired in for other stuff. Constraint is that the BAM stays in the aligner’s output order, i.e. name-grouped and not coordinate sorted.

The samtools equivalent is:

```
# ...
| samtools view -@ <threads> -b -o <sample>.<lane>.namegrouped.bam
```

Two of those flags are not the bwa defaults and TARS depends on both:

```text
-T 19   minimum alignment score to output. Lowered from the default 30 so short-anchor supplementaries at
splice junctions survive; TARS merges them back into the primary as junctions.

-h 75   XA cap, the maximum alternate loci listed per read. Alternates are carried as XA tags rather than
secondary records, and ISOFOX reads XA to identify multi-mapped fragments.

```

### 2. TARS

Source code at `/opt/repos/hmftools/tars` (branch `AUS418-tars-star-replacement`).

Command to run:

```bash
java -jar tars.jar \
-sample <sample> \
-input_bam <sample>.namegrouped.bam \
-ref_genome /path/to/Homo_sapiens_assembly38.alt.masked.with_rna_contigs.fasta \
-ref_genome_version V38 \
-contig_sidecar /path/to/ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv \
-rna_unmap_regions /path/to/rna/38/rna_excluded_regions.38.tsv \
-bamtool /path/to/samtools \
-output_dir /path/to/output/ \
-threads <threads>
```

Outputs, coordinate sorted and indexed:

```text
<sample>.tars.bam + .bai    lifted genomic BAM, input to REDUX
<sample>.tars.summary.tsv   counts summary of what liftback did

```

TARS needs the same combined FASTA that was used for alignment, since it must read the transcript contigs it is lifting from. It does not need an Ensembl cache; the exon and junction annotation it uses comes from the contig sidecar.

### 3. REDUX

Standard REDUX, no RNA-specific flags. Note the reference genome here is the plain genomic FASTA, not the transcriptome-augmented one.

### 4. ISOFOX

No CLI change. Input BAM is the REDUX output instead of the STAR output. ISOFOX defaults to the bwa/TARS multi-mapping model, so nothing needs to be passed for the new arm. There is a -star_aligner flag that restores the old STAR handling if we ever need to reproduce a pre-v3.1 run; it should not be set in v3.1.

## Resource files

Located at `/opt/repos/pipeline-resources/reference_genome_rna`. For GRCh38, these files are:

- `Homo_sapiens_assembly38.alt.masked.with_rna_contigs.fasta`:
  - plus transcript contigs, plus its full bwa-mem2 index alongside (.0123, .amb, .ann, .bwt.2bit.64, .pac, .fai).
  - Used by bwa-mem2 and by TARS.
- `ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv`:
  -contig sidecar, maps each transcript contig interval back to its genomic exon span. 
  - TARS only.
- `hmftools/dna/mappability/unmap_regions.dna.38.tsv`: 
  - TARS only.

GRCh37 files are also available at the same directory shape.

## Oncoanalyser changes

### Subworkflow

Subworkflows refer to those at `oncoanalyser/subworkflows/local`

The relevant subworkflows are READ_ALIGNMENT_DNA and READ_ALIGNMENT_RNA. 

Because bwa-mem2 will be used for both DNA and RNA alignment, workflow READ_ALIGNMENT_RNA will be removed and 
READ_ALIGNMENT_DNA will be adapted to accept RNA inputs. The workflow will be renamed to READ_ALIGNMENT.

### Process 

The BWAMEM2_ALIGN process is located at `oncoanalyser/modules/local/bwa-mem2/mem/main.nf`:

Process will be used for both DNA and RNA alignment, just taking different args. The `input` directive of
the process will now take `is_rna` (a boolean), which will decide the (different) args passed to bwa-mem2. the bwa-mem2 
process takes only DNA fastqs or only RNA fastqs. 

There's no need to distinguish the DNA vs RNA BAMs at the process level. The DNA vs RNA separation will be at the top 
level workflows (WGTS, TARGETED, etc) - see below.

### Workflows

Workflows (aka top level workflows) refer to those in `oncoanalyser/workflows` (e.g. WGTS, TARGETED, etc).

Workflows will distinguish between DNA and RNA alignment subworkflows via these imports: 

```nextflow
include { READ_ALIGNMENT as READ_ALIGNMENT_DNA }
include { READ_ALIGNMENT as READ_ALIGNMENT_RNA }
```
