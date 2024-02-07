FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/sigs/ && \
  wget -O /opt/sigs/sigs.jar 'https://github.com/hartwigmedical/hmftools/releases/download/sigs-v1.2/sigs_v1.2.jar'

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'openjdk >=8' && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
