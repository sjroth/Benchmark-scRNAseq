# Use miniconda3 as our base image.
FROM continuumio/miniconda3:4.10.3-alpine

# Configure Bioconda
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install latest version of STAR and clean-up miniconda.
RUN conda install -c bioconda star=2.7.10a -y && \
   conda clean --all --yes
