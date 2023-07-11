#!/usr/bin/env bash

python train_classifiers.py "data/test/inputs/3,7-dihydroxytropolone.gbk" "data/test/inputs/3,7-dihydroxytropolone.txt" \
    --output "data/test" \
    --no_SSN "True" \
    --antismash_version 5 \
    --rgi_version 5
