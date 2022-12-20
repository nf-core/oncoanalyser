FROM mambaorg/micromamba:0.24.0

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    chardet \
    cython \
    dask \
    mypy \
    numpy \
    'openjdk >=8' \
    pandas \
    scikit-allel \
    zlib && \
  micromamba clean --all --yes

USER root

RUN \
  apt-get update && \
  apt-get install -y build-essential git procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

RUN \
  pip install 'TransVar==2.5.8.20190813'

RUN \
  git clone https://github.com/hartwigmedical/peach.git -b v1.6 /opt/peach/

RUN \
  apt-get remove -y build-essential git && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

USER mambauser
