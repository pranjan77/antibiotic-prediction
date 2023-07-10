#!/usr/bin/env bash

INPUT_FILE=$1
INPUT_DIR=$PWD
OUTPUT_DIR=$2

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir ${OUTPUT_DIR}
fi

echo "Running antismash on $1"

antismash $INPUT_FILE \
    --fullhmmer \
    --genefinding-tool prodigal \
    --output-dir $OUTPUT_DIR
