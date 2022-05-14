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
# prerequisites
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y autoconf
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y pkg-config
RUN apt-get install -y libssl-dev
RUN apt-get install -y make
RUN apt-get install -y cmake
RUN apt-get install -y libhdf5-serial-dev

# install R
RUN wget https://cran.r-project.org/src/base/R-latest.tar.gz
RUN tar xvfz R-latest.tar.gz
RUN cd R-*/; ./configure; make; make install
# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron


