FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/teal/ && \
  wget -O /opt/teal/teal.jar 'https://github.com/hartwigmedical/hmftools/releases/download/teal-v1.0.1/teal-1.0.1.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    picard && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
