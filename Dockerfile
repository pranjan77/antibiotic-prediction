FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

ARG INS_COMMIT="94a47d68b9eb03425fc26690e0d3e2af13dd3524"

RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*


RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Setup mamba
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
RUN bash Mambaforge-Linux-x86_64.sh -b -p "/mambaforge"

# Download setup scripts
RUN wget https://raw.githubusercontent.com/dileep-kishore/antibiotic-prediction/$INS_COMMIT/setup/install_antismash.sh -P /deps
RUN wget https://raw.githubusercontent.com/dileep-kishore/antibiotic-prediction/$INS_COMMIT/setup/install_rgi.sh -P /deps
RUN wget https://raw.githubusercontent.com/dileep-kishore/antibiotic-prediction/$INS_COMMIT/setup/install_natural_product.sh -P /deps

# Run setup scripts
RUN bash /deps/install_antismash.sh
RUN bash /deps/install_rgi.sh $INS_COMMIT
RUN bash /deps/install_natural_product.sh $INS_COMMIT


# Clone antibiotic-prediction repo
ARG RUN_COMMIT="e29f7ebf97f24645f24077a70a34b92f01cc8cea"
RUN echo '12' >/dev/null && mkdir /deps && cd /deps && \
       git clone --branch main https://github.com/dileep-kishore/antibiotic-prediction.git && \
       cd antibiotic-prediction && git checkout $RUN_COMMIT

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
