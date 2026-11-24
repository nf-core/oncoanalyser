# Software versions

The main software used by this version of `oncoanalyser` is listed below. Version changes between releases are recorded
in the [CHANGELOG](https://github.com/nf-core/oncoanalyser/blob/master/CHANGELOG.md), and the exact versions used in
each run are written to `pipeline_info/nf_core_oncoanalyser_software_mqc_versions.yml` in the output directory.

## WiGiTS tools

Tools from [WiGiTS](https://github.com/hartwigmedical/hmftools), the open-source suite of genomic analysis tools developed
by the Hartwig Medical Foundation. These tools carry out most of the analysis in `oncoanalyser`.

| Tool               | Version | Description                                            |
| ------------------ | ------- | ------------------------------------------------------ |
| `AMBER`            | 4.3     | B-allele frequency analysis for purity and ploidy      |
| `BamTools`         | 1.6.1   | Alignment and coverage metrics                         |
| `CHORD`            | 2.1.2   | Homologous recombination deficiency (HRD) prediction   |
| `CIDER`            | 1.2     | Immune receptor (IG/TCR) CDR3 sequence detection       |
| `COBALT`           | 3.0     | GC-normalised read depth ratios                        |
| `CUPPA`            | 2.5.1   | Tissue of origin prediction                            |
| `ESVEE`            | 2.0.1   | Structural variant calling                             |
| `FastqTools`       | 1.0     | UMI extraction from FASTQ read names                   |
| `ISOFOX`           | 2.0.1   | RNA transcript expression and fusion analysis          |
| `LILAC`            | 2.0     | HLA class I typing and allele-specific status          |
| `LINX`             | 2.3.1   | SV clustering, fusions and driver event interpretation |
| `NEO`              | 1.3     | Neoepitope prediction                                  |
| `ORANGE`           | 5.0.1   | Summary report                                         |
| `PAVE`             | 1.9     | Gene impact annotation of small variants               |
| `PEACH`            | 2.0.0   | Pharmacogenomic haplotype calling                      |
| `PURPLE`           | 4.4     | Purity, ploidy and copy number estimation              |
| `QSEE`             | 1.0     | Sequencing and variant calling QC metrics              |
| `REDUX`            | 2.0.5   | Duplicate marking, UMI consensus and read unmapping    |
| `SAGE`             | 5.0.2   | SNV, MNV and INDEL calling                             |
| `SIGS`             | 1.2.1   | Mutational signature fitting                           |
| `TEAL`             | 1.4     | Telomere length and rearrangement characterisation     |
| `VirusInterpreter` | 1.7.2   | Oncoviral integration interpretation                   |
| `WISP`             | 1.3.1   | Tumor fraction estimation (`purity_estimate` mode)     |

## Other tools

External tools used for read alignment, read processing, oncoviral detection and reporting.

| Tool                            | Version | Description                              |
| ------------------------------- | ------- | ---------------------------------------- |
| `BWA-MEM2`                      | 2.3     | DNA read alignment                       |
| `fastp`                         | 0.23.4  | FASTQ splitting and UMI extraction       |
| `GATK4` (Picard MarkDuplicates) | 4.6.1.0 | RNA duplicate marking                    |
| `GRIDSS` (VIRUSbreakend)        | 2.13.2  | Oncoviral detection                      |
| `linxreport`                    | 1.2.0   | Interactive LINX report                  |
| `MultiQC`                       | 1.35    | Aggregated QC report                     |
| `RepeatMasker`                  | 4.1.5   | Repeat annotation for VIRUSbreakend      |
| `sambamba`                      | 1.0.1   | Conversion and sorting of DNA alignments |
| `samtools`                      | 1.21    | BAM sorting and indexing                 |
| `STAR`                          | 2.7.3a  | RNA read alignment                       |
