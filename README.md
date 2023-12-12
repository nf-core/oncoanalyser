# ![nf-core/oncoanalyser](docs/images/nf-core-oncoanalyser_logo_light.png#gh-light-mode-only) ![nf-core/oncoanalyser](docs/images/nf-core-oncoanalyser_logo_dark.png#gh-dark-mode-only)

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/oncoanalyser/results)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.5-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/launch%20on-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/oncoanalyser)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23oncoanalyser-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/oncoanalyser)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/oncoanalyser** is a Nextflow implementation of the comprehensive cancer DNA and RNA analysis and reporting
workflow from the Hartwig Medical Foundation. For detailed information on each component of the Hartwig Medical
Foundation workflow, please refer to [hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools/).

The oncoanalyser pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across
multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation
trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
implementation of this pipeline uses one container per process which makes it much easier to maintain and update
software dependencies. Where possible, these processes have been submitted to and installed from
[nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to
everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud
infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on
real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other
analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core
website](https://nf-co.re/oncoanalyser/results).

## Pipeline summary

The following processes and tools can be run with oncoanalyser:

* SNV and MNV calling (`SAGE`, `PAVE`)
* SV calling (`SV Prep`, `GRIDSS`, `GRIPSS`, `PURPLE`, `LINX`)
* CNV calling (`AMBER`, `COBALT`, `PURPLE`)
* Transcript analysis (`Isofox`)
* Oncoviral detection (`VIRUSBreakend`, `Virus Interpreter`)
* HLA calling (`LILAC`)
* HRD status prediction (`CHORD`)
* Mutational signature fitting (`Sigs`)
* Tissue of origin prediction (`CUPPA`)
* Report generation (`ORANGE`, `gpgr [LINX]`)

## Quick Start

Create a samplesheet containing your inputs:

```text
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
P1__wgts,P1,SA,tumor,dna,bam,/path/to/SA.tumor.dna.wgs.bam
P1__wgts,P1,SB,tumor,rna,bam,/path/to/SB.tumor.rna.wts.bam
P1__wgts,P1,SC,normal,dna,bam,/path/to/SC.normal.dna.wgs.bam
```

Launch oncoanalyser:

```bash
nextflow run umccr/oncoanalyser \
   -profile docker \
   --mode wgts \
   --genome GRCh38_hmf \
   --input samplesheet.csv \
   --outdir output/
```

## Documentation

The nf-core/oncoanalyser pipeline comes with documentation about the pipeline
[usage](https://nf-co.re/oncoanalyser/usage), [parameters](https://nf-co.re/oncoanalyser/parameters) and
[output](https://nf-co.re/oncoanalyser/output).

## Version support

As oncoanalyser is used in clinical settings and is subject to accreditation standards in some instances, there is a
need for long-term stability and reliability for feature releases in order to meet operational requirements. This is
accomplished through long-term support of several nominated feature releases, which all receive bug fixes and security
fixes during the period of extended support.

Each release that is given extended support is allocated a separate long-lived git branch with the 'stable' prefix, e.g.
`stable/1.2.x`, `stable/1.5.x`. Feature development otherwise occurs on the `main` branch.

## Credits

The oncoanalyser pipeline was written by Stephen Watts at the [University of Melbourne Centre for Cancer
Research](https://mdhs.unimelb.edu.au/centre-for-cancer-research).

We thank the following people for their extensive assistance in the development of this pipeline, listed in alphabetical
order:

* Charles Shale
* Oliver Hofmann
* Peter Priestley

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#oncoanalyser`
channel](https://nfcore.slack.com/channels/oncoanalyser) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

You can cite the oncoanalyser zenodo record for a specific version using the following doi:
[10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md)
file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia,
> Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
