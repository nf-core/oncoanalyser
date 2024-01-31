FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/amber/ && \
  wget -O /opt/amber/amber.jar 'https://github.com/hartwigmedical/hmftools/releases/download/amber-v4.0.rc1/amber_v4.0.rc1.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-copynumber \
    'openjdk >=8' \
    r-dplyr && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
