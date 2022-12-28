FROM docker.io/scwatts/gridss:2.13.2--3

ARG GH_BASE_URL='https://github.com/hartwigmedical/hmftools'
USER root

RUN \
  mkdir -p /opt/svprep/ && \
  wget -O /opt/svprep/svprep.jar 'https://github.com/hartwigmedical/hmftools/releases/download/sv-prep-v1.1/sv-prep_v1.1_beta.jar'

COPY assets/gridss.svprep /opt/svprep/gridss_svprep
RUN \
  chmod 755 /opt/svprep/gridss_svprep && \
  ln -s /opt/svprep/gridss_svprep /usr/local/bin/

USER mambauser
