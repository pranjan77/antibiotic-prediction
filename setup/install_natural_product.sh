#!/usr/bin/env bash

export PATH="/mambaforge/bin:$PATH"

# Set up natural product
mamba env create -f /deps/antibiotic-prediction/setup/env_natural_product.yml
