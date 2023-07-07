FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*


RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Clone directory
RUN git clone https://github.com/dileep-kishore/antibiotic-prediction.git

# Set up mambaforge
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
RUN bash Mambaforge-Linux-x86_64.sh -b

# FIXME: Install antismash 6? using conda
RUN mamba create -n antismash antismash

# Set up RGI5
RUN mamba env create -f antibiotic-prediction/env_rgi5.yml

# Set up natural product
RUN mamba env create -f antibiotic-prediction/env_natural_product.yml


#===========
