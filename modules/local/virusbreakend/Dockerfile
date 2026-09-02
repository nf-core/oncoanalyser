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
        'grep'

# Install Dfam database required for RepeatMasker
RUN \
    conda install -y curl && \
    curl -s https://www.dfam.org/releases/Dfam_3.7/families/Dfam_curatedonly.h5.gz | \
        gzip -cd > /env/share/RepeatMasker/Libraries/Dfam.h5

RUN \
    conda clean -yaf

# Move Conda environment into standard BioContainers base image
FROM quay.io/bioconda/base-glibc-busybox-bash:2.1.0

COPY --from=build /env/ /env/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
