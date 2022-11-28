FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/purple/ && \
  wget -O /opt/purple/purple.jar 'https://github.com/hartwigmedical/hmftools/releases/download/purple-v3.7/purple_v3.7.1.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bcftools \
    bioconductor-variantannotation \
    'circos >=0.69.6' \
    'openjdk >=8' \
    r-cairo \
    r-cowplot \
    r-dplyr \
    r-ggplot2 \
    r-tidyr \
    xorg-libxtst \
    zlib && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
