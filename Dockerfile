FROM ubuntu:16.04
MAINTAINER Chen Yuelong <yuelong.chen@aegicare.com>


ARG BUILD_PACKAGES="build-essential r-base r-base-core r-recommended r-api-3.5 libxml2-dev libcurl4-openssl-dev"
ARG DEBIAN_FRONTEND=noninteractive

# Update the repository sources list
RUN apt-get update && apt-get install --yes software-properties-common && \
    add-apt-repository --yes ppa:marutter/rrutter3.5 && \
    apt-get update && \
    apt-get install --yes \
              $BUILD_PACKAGES && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/BiocManager_1.30.4.tar.gz');BiocManager::install('DESeq2');BiocManager::install('pheatmap')" && \
    apt-get remove --purge --yes \
              build-essential git && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/*
ADD *.R /usr/local/bin/
RUN chmod +x /usr/local/bin/*

CMD R



