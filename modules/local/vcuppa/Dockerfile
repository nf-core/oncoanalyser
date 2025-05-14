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
        'hmftools-cuppa ==2.3.2' \
        'pytorch' \
        'torchvision' \
        'matplotlib-base'

RUN \
    conda clean -yaf

RUN \
    mkdir -p /opt/vcuppa/ && \
    wget -O /opt/vcuppa/vcuppa.jar https://github.com/hartwigmedical/hmftools/releases/download/v-cuppa-v0.1.0/v-cuppa_v0.1.0.jar

# Move Conda environment into standard BioContainers base image
FROM quay.io/bioconda/base-glibc-busybox-bash:3.1

COPY --from=build /env/ /env/
COPY --from=build /opt/vcuppa/ /opt/vcuppa/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
