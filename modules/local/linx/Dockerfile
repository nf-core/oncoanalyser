FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/linx/ && \
  wget -O /opt/linx/linx.jar 'https://github.com/hartwigmedical/hmftools/releases/download/linx-v1.24.1/linx_v1.24.1.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-gviz \
    'circos >=0.69.6' \
    'openjdk >=8' \
    'r-base == 4.1.3' \
    r-cowplot \
    r-dplyr \
    r-ggplot2 \
    r-magick \
    r-tidyr \
    xorg-libxtst && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
