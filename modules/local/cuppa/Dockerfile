FROM mambaorg/micromamba:0.24.0

USER root

RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/cuppa/ && \
  wget -O /opt/cuppa/cuppa.jar 'https://github.com/hartwigmedical/hmftools/releases/download/cuppa-v1.7/cuppa_v1.7.jar'

ARG URL_BASE=https://raw.githubusercontent.com/hartwigmedical/hmftools/cuppa-v1.7/cuppa/src/main/resources/cuppa-chart
ARG DIR_BASE=/opt/cuppa/chart
RUN \
  mkdir -p ${DIR_BASE}/ && \
    wget -P ${DIR_BASE}/ ${URL_BASE}/requirements.txt && \
    wget -P ${DIR_BASE}/ ${URL_BASE}/cuppa-chart.py && \
    wget -P ${DIR_BASE}/main/ ${URL_BASE}/main/create_chart.py && \
    wget -P ${DIR_BASE}/main/ ${URL_BASE}/main/prepare_data.py

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'matplotlib >=3.3,<4' \
    'numpy >=1.19,<2' \
    'openjdk >=8' \
    'pandas >=1.1,<2' \
    'scipy >=1,<2' \
    'seaborn >=0.11,<1' && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
