#!/usr/bin/env bash

GIT_COMMIT=$1

export PATH="/mambaforge/bin:$PATH"

# Set up natural product
wget https://raw.githubusercontent.com/dileep-kishore/antibiotic-prediction/$GIT_COMMIT/setup/env_natural_product.yml -P /deps
mamba env create -f /deps/env_natural_product.yml
