# Key concepts and conventions

## Alignment

BWA-MEM2 is used internally in `oncoanalyser` for alignment. The pipeline has been validated on and is compatible with
BAMs aligned with BWA-MEM, BWA-MEM2 and DRAGEN. Note that the mate CIGAR attribute is mandatory for any BAM records with
paired reads. Non-compatible BAMs may be rectified using tools such as the [Picard
FixMateInformation](https://gatk.broadinstitute.org/hc/en-us/articles/360036713471-FixMateInformation-Picard) routine.

## Unmapping problematic reads

After read alignment, [REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) is run and performs
'unmapping' of reads in pre-defined problematic regions which are discordant, have long soft clipping, or are in a
region of extreme high depth. The purpose of this unmapping step is to remove obvious poor alignments from the BAM prior
to running downstream tools. The unmapped reads are retained in the BAM. Overall, the problematic regions make up ~0.3%
of the genome and lead to ~3-6% of all reads being unmapped depending on genome version

## Deduplication, consensus and UMIs

In `oncoanalyser`, read deduplication is also performed by
[REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux). Deduplication aims to remove both PCR and
optical duplicates to avoid double counting of fragments. If UMIs (unique molecular identifiers) are present and
configured in `oncoanalyser`, then UMI aware deduplication will be performed. If duplicate fragments are found, then
REDUX marks all fragments as duplicates and creates a single consensus read with consensus bases and base qualities
computed. The consensus fragment is annotated as either single or dual strand. This allows downstream tools to
distinguish between high quality versus low quality consensus reads.

A detailed description of deduplication logic is available in the [REDUX
documentation](https://github.com/hartwigmedical/hmftools/tree/master/redux#deduplication).

## Error recalibration

Two types of sample specific error recalibration are currently performed in `oncoanalyser`. REDUX measures the rate of
microsatellite errors per consensus type (for UMIs: single vs dual stranded), repeat context, repeat length, and fits
these variables to a model (see [REDUX microsatellite jitter
modeling](https://github.com/hartwigmedical/hmftools/tree/master/redux#microsatellite-jitter-modelling) for details).
SAGE measures the rate of base errors per consensus type, trinucleotide context, and mutation type (see [SAGE
concepts](https://github.com/hartwigmedical/hmftools/tree/master/sage#key-concepts-in-sage) for details). In both cases
the rate of recalibrated errors are stored in lookup files and are used downstream in small variant calling.

## Gene and transcript definitions

Hartwig's universe of genes consists of all HGNC symbols with a matching Ensembl gene in GRCh38. The universe of
transcripts consists of all Ensembl transcripts belonging to a gene with a matching HGNC symbol. The canonical
transcript is set to the Ensembl canonical transcript. For GRCh37, the transcripts differ substantially as the Ensembl
database is no longer updated.

More details can be found on the [HMF Gene Utilities
documentation](https://github.com/hartwigmedical/hmftools/tree/master/gene-utils#overview-of-gene-configuration).

## Driver gene panel

The driver gene panel is a key configuration in `oncoanalyser`. Genes configured in this file are used to generate a BED
file which defines the `PANEL` tier for variant calling. Calling of driver events is controlled by the per gene
configuration in this file. Users may wish to modify the driver list to be more representative of their specific cancer
type (e.g. Adult driver genes are very different from pediatric).

Note that some reference data files in `oncoanalyser` have been generated based on the default pre-defined driver gene
panel. For example, driver likelihood estimates from
[PURPLE](https://github.com/hartwigmedical/hmftools/blob/master/purple/DriverCatalog.md#gene-driver-likelihood) were
trained on an adult pan-cancer cohort. Similarly, driver genes used in
[CUPPA](https://github.com/hartwigmedical/hmftools/tree/master/cuppa) were selected based on their presence/frequency in
this adult pan-cancer cohort. In the future we aim to make these reference data files more customisable.

More details on the driver gene panel can be found in the [PURPLE driver catalog
documentation](https://github.com/hartwigmedical/hmftools/blob/master/purple/DriverCatalog.md).
