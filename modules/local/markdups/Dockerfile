FROM docker.io/continuumio/miniconda3:23.10.0-1

RUN \
  conda install -y -n base conda-libmamba-solver && \
  conda config --set solver libmamba && \
  conda install -y -c bioconda -c conda-forge -c conda \
    'sambamba==1.0' 'samtools==1.17' 'openjdk >=8' && \
  conda clean -yaf

RUN \
  mkdir -p /opt/markdups/ && \
  wget -O /opt/markdups/markdups.jar 'https://github.com/hartwigmedical/hmftools/releases/download/mark-dups-v1.1/mark-dups_v1.1.rc1.jar'
