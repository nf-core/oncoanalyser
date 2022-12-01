FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y patch procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

ARG GH_BASE_URL=https://raw.githubusercontent.com/hartwigmedical/pipeline5
COPY assets/extractSigPredictHRD.R.patch /tmp/
RUN \
  mkdir -p /opt/chord/ && \
  wget -P /opt/chord/ ${GH_BASE_URL}/v5.30/cluster/images/extractSigPredictHRD.R && \
  patch -p1 /opt/chord/extractSigPredictHRD.R /tmp/extractSigPredictHRD.R.patch && \
  rm /tmp/extractSigPredictHRD.R.patch && \
  chmod 755 /opt/chord/extractSigPredictHRD.R

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-bsgenome \
    bioconductor-bsgenome.hsapiens.ucsc.hg19 \
    bioconductor-bsgenome.hsapiens.ucsc.hg38 \
    r-base \
    r-devtools \
    r-randomforest && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}:/opt/chord/"

RUN \
  R -e "devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor', ref='1.14')" && \
  R -e "devtools::install_github('https://github.com/UMCUGenetics/CHORD', ref='2.00')"
