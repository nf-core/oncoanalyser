FROM docker.io/continuumio/miniconda3:23.10.0-1

RUN \
  conda install -y -n base conda-libmamba-solver && \
  conda config --set solver libmamba && \
  conda install -y -c bioconda -c conda-forge -c conda \
    'bwa==0.7.17' \
    'sambamba==1.0' && \
  conda clean -yaf
