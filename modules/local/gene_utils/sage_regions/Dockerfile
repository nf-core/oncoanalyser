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
        'wget'

RUN \
    mkdir -p /opt/gene-utils/ && \
    wget -O /opt/gene-utils/gene-utils.jar 'https://github.com/hartwigmedical/hmftools/releases/download/gene-utils-v1.2/gene-utils_v1.2.jar'

RUN \
    conda clean -yaf

# Move Conda environment into standard BioContainers base image
FROM quay.io/bioconda/base-glibc-busybox-bash:3.1

COPY --from=build /env/ /env/
COPY --from=build /opt/gene-utils/ /opt/gene-utils/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
