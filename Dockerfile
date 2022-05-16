FROM rocker/tidyverse:4.2.0 as Rbase
RUN apt-get update && apt-get install -y libglpk-dev
# R packages that are not in conda
RUN R -e "install.packages('pacman', repos='http://cran.rstudio.com/')"
# SPONGE
RUN R -e "devtools::install_github('biomedbigdata/SPONGE')"
RUN R -e "if(!require(SPONGE)) stop('SPONGE not properly installed')"

# install psirc from git repository
ARG DEBIAN_FRONTEND=noninteractive
# system prerequisites
RUN apt-get update && apt-get install -y \
      apt-utils \
      autoconf \
      libcurl4-openssl-dev \
      pkg-config \
      libssl-dev \
      make \
      cmake \
      libhdf5-serial-dev \
      libbz2-dev \
      liblzma-dev \
      g++ \
      && rm -rf /var/lib/apt/lists/*

COPY R_p_install.R /
RUN Rscript /R_p_install.R \
      biomaRt \
      argparser \
      data.table \
      dplyr \
      ggplot2 \
      reshape2 \
      stringr \
      VennDiagram \
      Biostrings \
      MetBrewer \
      ensembldb \
      pheatmap \
      DESeq2 \
      EnhancedVolcano \
      doParallel \
      foreach \
      BSgenome \
      GenomicRanges \
      GenomicFeatures \
      seqinr

FROM nfcore/base:1.12.1
LABEL authors="Octavia Ciora, Leon Schwartz, Markus Hoffmann" \
      description="Docker image containing all software requirements for the nf-core/circrnasponging pipeline"
# add R and base software
COPY --from=Rbase . ./
# install psirc
COPY install_psirc.sh /
RUN bash /install_psirc.sh
# install PITA
COPY install_pita.sh /
RUN bash /install_pita.sh
# add psirc and PITA to PATH
ENV PATH /ext/psirc/psirc-quant/release/src:/ext/PITA:$PATH

# Install the conda environment
# All conda/bioconda dependencies are listed there
COPY environment.yml /

RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-circrnasponging/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-circrnasponging > nf-core-circrnasponging.yml
# include firefox and geckodriver
ENV GECKODRIVER_VER v0.29.0
ENV FIREFOX_VER 87.0
# install dependencies
RUN set -x \
      && echo "deb http://deb.debian.org/debian/ unstable main contrib non-free" >> /etc/apt/sources.list.d/debian.list \
      && apt-get update \
      && apt-get upgrade -y \
      && apt-get install -y \
            firefox-esr

# Add latest FireFox
RUN set -x \
      && apt-get install -y \
            ibx11-xcb1 \
            libdbus-glib-1-2 \
      && curl -sSLO https://download-installer.cdn.mozilla.net/pub/firefox/releases/${FIREFOX_VER}/linux-x86_64/en-US/firefox-${FIREFOX_VER}.tar.bz2 \
      && tar -jxf firefox-* \
      && mv firefox /opt/ \
      && chmod 755 /opt/firefox \
      && chmod 755 /opt/firefox/firefox
  
# Add geckodriver
RUN set -x \
      && curl -sSLO https://github.com/mozilla/geckodriver/releases/download/${GECKODRIVER_VER}/geckodriver-${GECKODRIVER_VER}-linux64.tar.gz \
      && tar zxf geckodriver-*.tar.gz \
      && mv geckodriver /usr/bin/