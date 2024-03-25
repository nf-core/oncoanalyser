FROM continuumio/miniconda3:23.5.2-0-alpine as build

RUN \
  conda install conda-libmamba-solver

RUN \
  echo -e > ~/.condarc '\
solver: libmamba\n\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
  - defaults'

RUN \
  conda create -y -p /env/ \
    'gridss=2.13.2=h50ea8bc_3' \
    'bash' \
    'coreutils' \
    'findutils' \
    'gawk' \
    'grep' \
    'procps-ng' \
    'time' \
    'util-linux' \
    'which'

RUN \
  conda clean -yaf

# Move Conda environment into distroless image
FROM gcr.io/distroless/base-debian11:latest

COPY --from=build /env/ /env/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"

# Symlink system executables as required by VIRUSBreaked
SHELL ["/env/bin/bash", "-c"]
RUN \
  ln -s /env/bin/bash /bin/bash && \
  ln -s /env/bin/env /usr/bin/env && \
  ln -s /env/bin/time /usr/bin/time
