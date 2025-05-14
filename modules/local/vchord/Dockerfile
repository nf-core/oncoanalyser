FROM continuumio/miniconda3:25.3.1-1 as build

RUN \
    echo -e > ~/.condarc '\
solver: libmamba\n\
channels:\n\
    - conda-forge\n\
    - bioconda\n\
    - defaults'

RUN \
    conda create -y -p /env/ \
        'openjdk >=8,<=21' \
        'matplotlib-base' \
        'numpy' \
        'pandas' \
        'pytorch' \
        'torchvision'

RUN \
    conda clean -yaf

RUN \
    mkdir -p /opt/vchord/ && \
    wget -O /opt/vchord/vchord.jar https://github.com/hartwigmedical/hmftools/releases/download/v-chord-v1.0/v-chord_v1.0.jar

# Move Conda environment into standard BioContainers base image
FROM quay.io/bioconda/base-glibc-busybox-bash:3.1

COPY --from=build /env/ /env/
COPY --from=build /opt/vchord/ /opt/vchord/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
