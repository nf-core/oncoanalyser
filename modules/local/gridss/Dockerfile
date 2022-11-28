FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'gridss==2.13.2' && \
  micromamba clean --all --yes

RUN \
  sed -i \
    -e '209s/grep samtools/head -n1 | grep samtools/' \
    /opt/conda/share/gridss-2.13.2-1/gridss_extract_overlapping_fragments

RUN \
  ln -sf \
    /opt/conda/pkgs/gridss-2.13.2-h20b1175_1/share/gridss-2.13.2-1/gridss \
    /opt/conda/bin/gridss

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
