FROM docker.io/scwatts/gridss:2.13.2--3

ARG GH_BASE_URL='https://github.com/hartwigmedical/hmftools'
USER root

RUN \
  apt-get update && \
  apt-get install -y patch && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN \
  mkdir -p /opt/svprep/ && \
  wget -O /opt/svprep/svprep.jar 'https://github.com/hartwigmedical/hmftools/releases/download/sv-prep-v1.2.3/sv-prep_v1.2.3.jar' && \
  wget -O /opt/svprep/gridss_svprep 'https://raw.githubusercontent.com/hartwigmedical/hmftools/sv-prep-v1.2.3/sv-prep/src/main/resources/gridss/gridss.run.sh'

COPY assets/gridss_svprep.patch /tmp/
RUN \
  patch -lp1 -i /tmp/gridss_svprep.patch /opt/svprep/gridss_svprep && \
  chmod 755 /opt/svprep/gridss_svprep && \
  ln -s /opt/svprep/gridss_svprep /usr/local/bin/

USER mambauser
