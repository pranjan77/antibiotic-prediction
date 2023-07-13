#!/usr/bin/env bash

mkdir -p /tmp
cd /tmp

apt-get update && \
    apt-get -y install apt-transport-https gnupg

wget http://dl.secondarymetabolites.org/antismash-stretch.list -O /etc/apt/sources.list.d/antismash.list && \
    wget -q -O- http://dl.secondarymetabolites.org/antismash.asc | apt-key add -

apt-get update && \
    apt-get -y install hmmer2 hmmer diamond-aligner fasttree prodigal ncbi-blast+ muscle glimmerhmm

export PATH="/mambaforge/bin:$PATH"

mamba create -n antismash python=3.9
source activate antismash

# Alternative to line 12
# mamba install -c bioconda hmmer2 hmmer diamond=2.0.14 fasttree prodigal blast muscle glimmerhmm

wget https://dl.secondarymetabolites.org/releases/7.0.0/antismash-7.0.0.tar.gz  \
      && tar -zxf antismash-7.0.0.tar.gz

pip install ./antismash-7.0.0
download-antismash-databases && antismash --check-prereqs

# pip3 install nose jinja2
# pip3 install jsonrpcbase numpy pandas biopython
# pip3 install coverage
