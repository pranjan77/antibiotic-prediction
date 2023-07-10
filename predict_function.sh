#!/usr/bin/env bash

GENOME=$1
GENOME_NAME=$(basename "${GENOME%%.*}")
OUTPUT_DIR=$2/$GENOME_NAME
mkdir -p $OUTPUT_DIR

# Step 1: Run antismash
echo "Running antismash"
if [ ! -d ${OUTPUT_DIR/antismash} ]; then
    mkdir -p $OUTPUT_DIR/antismash
    bash run_antismash5.sh $GENOME $OUTPUT_DIR/antismash
fi

# Step 2: For each BGC get fasta and run RGI
source deactivate
source activate rgi5
echo "Running RGI"
if [ ! -d ${OUTPUT_DIR/rgi} ]; then
    mkdir -p $OUTPUT_DIR/rgi
    for f in $OUTPUT_DIR/antismash/*region*.gbk; do
        BGC=$(basename "$f" .gbk)
        mkdir $OUTPUT_DIR/rgi/$BGC
        python gbk2fasta.py "$f" $OUTPUT_DIR/rgi/"$BGC".fna
        rgi main -i $OUTPUT_DIR/rgi/"$BGC".fna -o $OUTPUT_DIR/rgi/"$BGC"/"$BGC" --include_loose
    done
fi

# Step 3: Run prediction script for each BGC
source deactivate
source activate natural_product
echo "Running BGC activity prediction"
python predict_function.py $OUTPUT_DIR/antismash $OUTPUT_DIR/rgi \
    --data_dir data \
    --output_dir $OUTPUT_DIR \
    --classifiers tree \
    --antismash_version 5 \
    --rgi_version 5
