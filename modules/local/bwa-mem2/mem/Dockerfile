FROM docker.io/mambaorg/micromamba:1.5.6

USER root

RUN \
    apt-get update && \
    apt-get install -y procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER mambauser

RUN \
    micromamba install -y -n base -c bioconda -c conda-forge \
        'bwa-mem2 ==2.2.1' \
        'sambamba ==1.0' && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
