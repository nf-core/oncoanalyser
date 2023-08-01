FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/cobalt/ && \
  wget -O /opt/cobalt/cobalt.jar 'https://github.com/hartwigmedical/hmftools/releases/download/cobalt-v1.15.2/cobalt_v1.15.2.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-copynumber \
    'openjdk >=8' \
    r-dplyr && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
