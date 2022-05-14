FROM r-base:4.2.0 as Rbase
# Add R and instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
# R packages that are not in conda
RUN R -e "install.packages(c('pacman'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "pacman::p_load(SPONGE, biomaRt, argparser, data.table, dplyr, ggplot2, reshape2, stringr, VennDiagram, Biostrings, MetBrewer, ensembldb, pheatmap, DESeq2, EnhancedVolcano, doParallel, foreach, BSgenome, GenomicRanges, GenomicFeatures, seqinr)"


FROM nfcore/base:1.12.1
LABEL authors="Octavia Ciora, Leon Schwartz, Markus Hoffmann" \
      description="Docker image containing all software requirements for the nf-core/circrnasponging pipeline"

# install psirc from git repository
ARG DEBIAN_FRONTEND=noninteractive
# psirc prerequisites
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
      g++
      # && rm -rf /var/lib/apt/lists/*
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

# add R
COPY --from=Rbase . ./
