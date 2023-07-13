#!/usr/bin/env bash

GIT_COMMIT=$1

export PATH="/mambaforge/bin:$PATH"

# Set up RGI5
wget https://raw.githubusercontent.com/dileep-kishore/antibiotic-prediction/$GIT_COMMIT/setup/env_rgi5.yml -P /deps
mamba env create -f /deps/env_rgi5.yml
