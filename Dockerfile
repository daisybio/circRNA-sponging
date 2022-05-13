FROM nfcore/base:1.12.1
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

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

# install R
RUN apt-get update && apt-get install -y r-base

# R packages that are not in conda
RUN R -e "install.packages(c('pacman'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "pacman::p_load(SPONGE, biomaRt, argparser, data.table, dplyr, ggplot2, reshape2, stringr, VennDiagram, Biostrings, MetBrewer, ensembldb, pheatmap, DESeq2, EnhancedVolcano, doParallel, foreach, BSgenome, GenomicRanges, GenomicFeatures, seqinr)"

# install psirc from git repository
ARG DEBIAN_FRONTEND=noninteractive
# psirc prerequisites
RUN apt-get install -y apt-utils
RUN apt-get install -y autoconf
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y pkg-config
RUN apt-get install -y libssl-dev
RUN apt-get install -y make

RUN apt-get install -y cmake
RUN apt-get install -y libhdf5-serial-dev
# install psirc
COPY install_psirc.sh /
RUN bash /install_psirc.sh
# install PITA
COPY install_pita.sh /
RUN bash /install_pita.sh

COPY . ./