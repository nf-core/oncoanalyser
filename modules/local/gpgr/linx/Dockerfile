FROM docker.io/scwatts/gpgr:1.4.3

RUN \
  mamba install -y patch && mamba clean --all --force-pkgs-dirs --yes

COPY assets/linx.R.patch /tmp/

RUN \
  git clone https://github.com/umccr/gpgr.git -b v1.4.4 /tmp/gpgr/ && \
  patch /tmp/gpgr/R/linx.R /tmp/linx.R.patch && \
  R CMD INSTALL --build /tmp/gpgr && \
  rm -rf /tmp/{gpgr,gpgr_1.4.4_R_x86_64-conda-linux-gnu.tar.gz,linx.R.patch}
