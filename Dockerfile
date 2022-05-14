FROM nfcore/base:1.12.1 AS nfcore
LABEL authors="Octavia Ciora, Leon Schwartz, Markus Hoffmann" \
      description="Docker image containing all software requirements for the nf-core/circrnasponging pipeline"

# Install the conda environment
# All conda/bioconda dependencies are listed there
COPY environment.yml /

RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-circrnasponging/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-circrnasponging > nf-core-circrnasponging.yml


ARG DEBIAN_FRONTEND=noninteractive
# prerequisites for psirc and PITA
RUN apt-get update -qq && \
      apt-get install -y apt-utils \
      autoconf \
      libcurl4-openssl-dev \
      pkg-config \
      libssl-dev \
      make \
      cmake \
      libhdf5-serial-dev \
      # Install system dependencies for R
      apt-transport-https \
      build-essential \
      curl \
      gfortran \
      libatlas-base-dev \
      libbz2-dev \
      libcairo2 \
      libicu-dev \
      liblzma-dev \
      libpango-1.0-0 \
      libpangocairo-1.0-0 \
      libreadline-dev \
      libpcre3-dev \
      libtcl8.6 \
      libtiff5 \
      libtk8.6 \
      libx11-6 \
      libxt6 \
      locales \
      tzdata \
      zlib1g-dev
    
# Install system dependencies for the tidyverse R packages
RUN apt-get install -y \
      libssl-dev \
      pandoc \
      libxml2-dev

# install R
RUN wget https://cran.r-project.org/src/base/R-latest.tar.gz
RUN tar xvfz R-latest.tar.gz
RUN cd R-*/; ./configure; make; make install
# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron


