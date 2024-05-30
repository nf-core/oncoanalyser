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
        'openjdk >=8' \
        'sambamba >=1.0.1' \
        'samtools >=1.17'

RUN \
    conda clean -yaf

RUN \
    wget -P /tmp/ https://github.com/hartwigmedical/hmftools/releases/download/mark-dups-v1.1.6_beta/mark-dups_v1.1.6_beta.jar

# Move Conda environment into standard BioContainers base image
FROM quay.io/bioconda/base-glibc-busybox-bash:2.1.0

COPY --from=build /env/ /env/
COPY --from=build /tmp/mark-dups_v1.1.6_beta.jar /opt/markdups/markdups.jar

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
