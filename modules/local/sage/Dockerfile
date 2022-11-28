FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/sage/ && \
  wget -O /opt/sage/sage.jar 'https://github.com/hartwigmedical/hmftools/releases/download/sage-v3.2/sage_v3.2.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-genomicranges \
    bioconductor-variantannotation \
    font-ttf-dejavu \
    'openjdk >=8' \
    r-dplyr \
    r-ggplot2 \
    r-tidyr \
    xorg-libxt \
    zlib && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

