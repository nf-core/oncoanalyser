FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/cuppa/ && \
  wget -O /opt/cuppa/cuppa.jar 'https://github.com/hartwigmedical/hmftools/releases/download/cuppa-v1.8.1/cuppa_v1.8.1.jar'

ARG URL_BASE=https://raw.githubusercontent.com/hartwigmedical/hmftools/cuppa-v1.8.1
ARG URL_BASE_CHART=${URL_BASE}/cuppa/src/main/resources/cuppa-chart
ARG DIR_BASE=/opt/cuppa
RUN \
  mkdir -p ${DIR_BASE}/{chart,report} && \
    wget -P ${DIR_BASE}/chart/ ${URL_BASE_CHART}/requirements.txt && \
    wget -P ${DIR_BASE}/chart/ ${URL_BASE_CHART}/cuppa-chart.py && \
    wget -P ${DIR_BASE}/chart/main/ ${URL_BASE_CHART}/main/create_chart.py && \
    wget -P ${DIR_BASE}/chart/main/ ${URL_BASE_CHART}/main/prepare_data.py

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'matplotlib >=3.3,<4' \
    'numpy >=1.19,<2' \
    'openjdk >=8' \
    'pandas >=1.1,<2' \
    'scipy >=1,<2' \
    'seaborn >=0.11,<1' \
    'r-base ==4.2.2' \
    'r-cowplot' \
    'r-ggplot2' \
    'r-gridextra' \
    'r-stringi' \
    'r-tidyr' && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
