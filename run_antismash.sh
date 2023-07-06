#!/usr/bin/env bash

INPUT_FILE=$1
INPUT_DIR=$PWD
OUTPUT_DIR=$2

# Links within the container
readonly CONTAINER_SRC_DIR=/input
readonly CONTAINER_DST_DIR=/output

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir ${OUTPUT_DIR}
fi

echo $1

sudo docker run \
    --volume "$PWD:/input:ro" \
    --volume "$PWD:/output:rw" \
    --detach=false \
    --rm \
    --user=$(id -u):$(id -g) \
    antismash/standalone:7.0.0 \
    $INPUT_FILE \
    --fullhmmer \
    --genefinding-tool prodigal \
    --output-dir $OUTPUT_DIR
