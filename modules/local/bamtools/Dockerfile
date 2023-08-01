FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/bamtools/ && \
  wget -O /opt/bamtools/bamtools.jar 'https://github.com/hartwigmedical/hmftools/releases/download/bam-tools-v1.1/bam-tools_v1.1.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'openjdk >=8' && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
