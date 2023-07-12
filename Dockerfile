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

RUN echo '12' >/dev/null && mkdir /deps && cd /deps && \
       git clone --branch main https://github.com/dileep-kishore/antibiotic-prediction.git

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
RUN bash Mambaforge-Linux-x86_64.sh -b -p "/mambaforge"

RUN bash /deps/antibiotic-prediction/setup/install_antismash.sh
RUN bash /deps/antibiotic-prediction/setup/install_rgi.sh
RUN bash /deps/antibiotic-prediction/setup/install_natural_product.sh

# -----------------------------------------

#ENV PATH=$PATH1

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
