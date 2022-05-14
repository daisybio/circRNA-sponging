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

# install R
ARG R_VERSION=4.2.0
RUN wget https://cdn.rstudio.com/r/ubuntu-1604/pkgs/r-${R_VERSION}_1_amd64.deb
RUN apt-get install -y gdebi-core
RUN gdebi r-${R_VERSION}_1_amd64.deb
# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron


